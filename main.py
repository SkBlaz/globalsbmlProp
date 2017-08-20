## this is the main analysis class

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
         # Print function arguments. 
 #         if (math.getNumChildren() > 1):
 #             print("(" + (math.getLeftChild()).getName());
 # # Print function arguments. 

 #             for n in range (1, math.getNumChildren()):
 #                 try:
 #                     print(", " + (math.getChild(n)).getName());
 #                 except:
 #                     pass
 
 #         print(") := ");
 
         # # Print function body. 
         if (math.getNumChildren() == 0):
             print("(no body defined)");
         else:
             math = math.getChild(math.getNumChildren() - 1);
             formula = formulaToString(math);
             return formula
#             print(formula + "\n");

def getModelMath(genModels,cmprt='all'):

    ## this function obtaines reaction and other math related to specific models and saves them into a dataframe
    
    cpfor = defaultdict(list)
    for candidate in genModels:
        document = readSBML(candidate);
        if (document.getNumErrors() > 0):
            pass
        else:
            model = document.getModel();

            formulas = []
            for n in range(0,model.getNumFunctionDefinitions()):
                formulas.append(printFunctionDefinition(n + 1, model.getFunctionDefinition(n)))
                
            for i in range(0, model.getNumCompartments()):
                sp = model.getCompartment(i)
                if cmprt != 'all':
                    if sp.getId() in cmprt:
                        cpfor[sp.getId()].append(formulas)
                else:
                    cpfor[sp.getId()].append(formulas)

    ## get compartment|formula structure
    return cpfor


def intra_compartment_distances(formula_file):

    import editdistance as ed
    from itertools import combinations
#    editdistance.eval('banana', 'bahama')
    distobj = {}
    
    for k,v in formula_file.items():

        all_formulas = []
        total_combinations=0
        current_ed = 0
        for flist in v:
            for formula in flist:
                all_formulas.append(formula)

        for f1,f2 in combinations(all_formulas,2):
            current_ed += ed.eval(f1,f2)
            total_combinations += 1

        try:
            distobj[k] = float(current_ed/total_combinations)
        except:
            pass

    ## perhaps plot this here..
    print(distobj)

def inter_component_distances(formula_file):

    import editdistance as ed
    from itertools import combinations
    import numpy as np
    import seaborn as sns
    
    print("Starting pairwise distance measurements..")
    distframe = pd.DataFrame()
    ## double loop for pairwise distances v is of form list of lists
    partial = 0
    totlen = len(formula_file.keys())
    for k,v in formula_file.items():
        partial+=1
        if partial % 10 == 0:
            print(float(partial*100/totlen),"%","complete.")
        for k2,v2 in formula_file.items():
            
            ## first get representative formulas for individual components
            v1_rep = [formula for sublist in v for formula in sublist]
            v2_rep = [formula for sublist in v2 for formula in sublist]                    
            distMinAvg = np.mean([ed.eval(s1,s2) for s1 in v1_rep for s2 in v2_rep])

            if distMinAvg >= 0:
                distframe = distframe.append({'component1' : k, 'component2' : k2, 'distance' : distMinAvg},ignore_index=True)

    indata = distframe.pivot("component1","component2","distance")
    ax = sns.heatmap(indata)
    plt.show()

if __name__ == "__main__":

    import argparse    
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats",help="Some basic statistics")
    parser.add_argument("--math",help="Math based process clustering")
    parser.add_argument("--lev",help="Distances within individual components")
    parser.add_argument("--interlev",help="Distances within individual components")       
    args = parser.parse_args()

    print("Beginning extraction..")
    
    datafolder = "data/BioModels_Database-r31_pub-sbml_files/curated"
    model_getter = model_generator(datafolder)    

    ## query parts of the cell/ organism
    compartments_to_check=['cell','nucleus','plasma','nuclei','CellSurface','cytosol','vacuole','Lysosome','Mitochondria','cellsurface','Endosome']
    
    if args.stats:
        ## those are some basic numeric statistics regarding individual models
        get_basic_stats(model_getter,compartment=compartments_to_check)

    comp_formulas = getModelMath(model_getter,cmprt='all')

    ## this computes inter-component distances, pairwise. even within a single component, there is large variability present within individual component!
    if args.interlev:
        ## this only gets the saved data, which is further
        inter_component_distances(comp_formulas)
