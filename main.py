
## this is the main analysis class

import re
import multiprocessing as mp
import itertools
import xml.etree.ElementTree as ET
import sys
import os.path
from libsbml import *  ## the main library
import glob ## for file search
from collections import Counter,defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def printAnnotation(sb, id=""):
    if (not sb.isSetAnnotation()):
        return;        
   
    pid = "";
    
    if (sb.isSetId()):
        pid = sb.getId();
    print("----- " + sb.getElementName() + " (" + pid
                     + ") annotation -----" + "\n");
    print(sb.getAnnotationString() + "\n");
    print("\n");


def model_generator(infolder):
    candidate_models = glob.glob(infolder+"/*.xml") ## get only xml files
    for cand in candidate_models:
        yield cand

def get_datapoint(model,sp):
    return {"FunctionDefinitions" : model.getNumFunctionDefinitions(),
            "unitDefinitions" : model.getNumUnitDefinitions(),
            "speciesTypes" : model.getNumSpeciesTypes(),
            "compartments" : model.getNumCompartments(),
            "species" : model.getNumSpecies(),
            "parameters" : model.getNumParameters(),
            "initialAssignments" : model.getNumInitialAssignments(),
            "rules" : model.getNumRules(),
            "constraints" : model.getNumConstraints(),
            "reactions" : model.getNumReactions(),
            "events" : model.getNumEvents(),"Compartment" : sp.getId()}

def get_basic_stats(sbmlfolder,compartment="all"):

    ## read and print the basic stats for individual .xml entries

#    candidate_models = glob.glob(sbmlfolder+"/*.xml") ## get only xml files
    all_levels = []
    df = pd.DataFrame()  ## this is the container for the basic stats
    all_compartments = []
    for candidate in sbmlfolder:
        document = readSBML(candidate);
        if (document.getNumErrors() > 0):
            pass
        else:
            level = document.getLevel();
            all_levels.append(level)
            model = document.getModel();
            if model !=None:
                for i in range(0, model.getNumCompartments()):
                    sp = model.getCompartment(i);
                    if compartment != "all":
                        ## add data only if specific compartment
                        all_compartments.append(sp.getId())
                        if sp.getId() in compartment:
                            datapoint = get_datapoint(model,sp)
                            df = df.append(datapoint,ignore_index=True)    
                    else:
                        datapoint = get_datapoint(model,sp)
                        df = df.append(datapoint,ignore_index=True)
                        
    ## analyze the obtained dataset
    print(Counter(all_levels))
#    print(set(all_compartments))
    print(df.describe()) ## this is to be further plotted
    g = sns.pairplot(df,hue="Compartment",vars=['reactions','species','constraints','FunctionDefinitions','rules'])
    plt.show()


def printReactionMath(n, r):
    if (r.isSetKineticLaw()):
        kl = r.getKineticLaw();
        if (kl.isSetMath()):
            formula = formulaToString(kl.getMath());
            print("Reaction " + str(n) + ", formula: " + formula + "\n");

def printFunctionDefinition(n, fd):
     if (fd.isSetMath()):
 #        print("FunctionDefinition " + str(n) + ", " + fd.getId());
 
         math = fd.getMath();
         # # Print function body. 
         if (math.getNumChildren() == 0):
             print("(no body defined)");
         else:
             math = math.getChild(math.getNumChildren() - 1);
             formula = formulaToString(math);
             return formula
#             print(formula + "\n");

def getModelMath(genModels,cmprt="all",goterms="all"):

    ## this function obtaines reaction and other math related to specific models and saves them into a dataframe
    
    cpfor = defaultdict(list)
    go_dict = defaultdict(list)
    for candidate in genModels:
        document = readSBML(candidate);
        if (document.getNumErrors() > 0):
            pass
        else:
            model = document.getModel();

            try:
                psx = RDFAnnotationParser()
                terms = psx.parseCVTerms(model)
                annotation = terms.toXMLString()
                tree = ET.ElementTree(ET.fromstring(annotation))
                goterms = []
                for elem in tree.iter():
                    if "GO:" in str(elem.attrib):
                        m = re.search("GO:\d{7}",str(elem.attrib))
                        goterms.append(m.group(0))
                                               
            except:
                pass
            
            
            formulas = []
            for n in range(0,model.getNumFunctionDefinitions()):
                formulas.append(printFunctionDefinition(n + 1, model.getFunctionDefinition(n)))
            for el in goterms:
                for f in formulas:
                    go_dict[el].append(f)
                
            for i in range(0, model.getNumCompartments()):
                sp = model.getCompartment(i)
                if cmprt != 'all':
                    if sp.getId() in cmprt:
                        cpfor[sp.getId()].append(formulas)
                else:
                    cpfor[sp.getId()].append(formulas)

    ## get compartment|formula structure
    return (cpfor,go_dict)

def inter_component_distances(formula_file,measure="ED",precomputed=None):

    import seaborn as sns
    if precomputed == None:
        if measure == "ED":
            import editdistance as ed
        elif measure == "fuzzy" or measure == "fuzzy_plain":
            from fuzzywuzzy import fuzz
            from fuzzywuzzy import process
        else:
            pass

        import numpy as np

        pool = mp.Pool(mp.cpu_count())
        print("Starting pairwise distance measurements..")
        distframe = pd.DataFrame()
        ## double loop for pairwise distances v is of form list of lists
        partial = 0
        totlen = len(formula_file.keys())
        for k,v in formula_file.items():
            partial+=1
            if partial % 1 == 0:
                print(float(partial*100/totlen),"%","complete.")
            for k2,v2 in formula_file.items():
            
                ## first get representative formulas for individual components
                v1_rep = [formula for sublist in v for formula in sublist]
                v2_rep = [formula for sublist in v2 for formula in sublist]

                all_pairs = []
                for f1 in v1_rep:
                    for f2 in v2_rep:
                        all_pairs.append((f1,f2))
                
                if measure == "ED":
                    distMinAvg = np.mean([pool.apply(ed.eval, args=(x,y,)) for x,y in all_pairs])
                if measure == "fuzzy":
                    distMinAvg = np.mean([pool.apply(fuzz.partial_ratio, args=(x,y,)) for x,y in all_pairs])
                if measure == "fuzzy_plain":
                    distMinAvg = np.mean([pool.apply(fuzz.ratio, args=(x,y,)) for x,y in all_pairs])

                if distMinAvg >= 0:
                    distframe = distframe.append({'First component' : k, 'Second component' : k2, 'distance' : distMinAvg},ignore_index=True)
    else:
        distframe = pd.read_csv(precomputed)
                    
    indata = distframe.pivot("First component","Second component","distance")
    ax = sns.heatmap(indata,cmap="BuGn")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.show()

    tmpframe = distframe[distframe['First component'] == distframe['Second component']]
    tmpframe = tmpframe.sort_values(['distance'])
    sns.barplot(x="First component",y="distance",data=tmpframe)
    plt.xticks(rotation=90)
    plt.ylabel("Intra-compartment distance")
    plt.show()

    distframe.to_csv("distances"+measure+".csv")


def get_similarity_list(fname):
    distframe = pd.read_csv(fname)
    distframe= distframe[distframe['First component'] != distframe['Second component']].sort_values(['distance'],ascending=False)
    distframe.to_csv(fname.split(".")[0]+"similar_list.csv")
    distframe.to_latex(fname.split(".")[0]+"similar_list.tex")
        
    
if __name__ == "__main__":

    import argparse    
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats",help="Some basic statistics")
    parser.add_argument("--lev",help="Distances within individual components")
    parser.add_argument("--interlev",help="Distances within individual components")
    parser.add_argument("--interfuzzy",help="Distances within individual components")
    parser.add_argument("--simstats",help="Most similar, yet not the same")
    parser.add_argument("--interfuzzybasic",help="Basic fuzzy algorithm")
    parser.add_argument("--goterms",help="Use GO terms?")  
    args = parser.parse_args()

    print("Beginning extraction..")
    
    datafolder = "data/BioModels_Database-r31_pub-sbml_files/curated"
    model_getter = model_generator(datafolder)

    ## query parts of the cell/ organism
    compartments_to_check=['cell','nucleus','plasma','nuclei','CellSurface','cytosol','vacuole','Lysosome','Mitochondria','cellsurface','Endosome']
    
    if args.stats:
        ## those are some basic numeric statistics regarding individual models
        get_basic_stats(model_getter,compartment=compartments_to_check)

    compartment_formulas, go_formulas = getModelMath(model_getter,cmprt=compartments_to_check)

    if args.goterms:
        print(len(go_formulas.keys())," Individual GO terms found.")    
        comp_formulas = go_formulas
    else:
        comp_formulas = compartment_formulas
    
    if args.interlev:

        ## this only gets the saved data, which is further        
        inter_component_distances(comp_formulas)

    if args.interfuzzy:
        inter_component_distances(comp_formulas,measure="fuzzy")

    if args.interfuzzybasic:
        inter_component_distances(comp_formulas,measure="fuzzy_plain")

    if args.simstats:
        get_similarity_list(args.simstats)


