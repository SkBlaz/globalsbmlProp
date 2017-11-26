
## this is the main analysis class

import re
import numpy as np
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

pool = mp.Pool(mp.cpu_count())

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


def navigateMath(node):
  #
  # first take the current element
  #
  if node.getType() == AST_PLUS:
    if node.getNumChildren() == 0:
       print ("0")
       return
    navigateMath(node.getChild(0))
    for i in range(1, node.getNumChildren()):
      print ("+")
      navigateMath(node.getChild(i))
    
  elif node.getType() == AST_REAL:
      # this will be constants 
      print (node.getReal())
  elif node.getType() == AST_NAME:
      print (node.getName())
  else:
      # handle more cases here ...
      pass

def getModelMath(genModels,cmprt="all",goterms="all"):

    ## this function obtaines reaction and other math related to specific models and saves them into a dataframe
    
    cpfor = defaultdict(list)
    go_dict = defaultdict(list)
    go_dict_ast = defaultdict(list)
    feature_list = ["<apply>","<times>","<power>","<divide>","<cn>","<ci>","<plus>","<minus>"]

    for x in range(2,3,1):
        extended = ["".join(x) for x in itertools.permutations(feature_list, x)]
        feature_list+=extended
        
    print("Possibly using:",len(feature_list),"fingerprints..")

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
                form = printFunctionDefinition(n + 1, model.getFunctionDefinition(n))
                formulas.append(form)                

            
            for el in goterms:
                for f in formulas:

                    treeRep = writeMathMLToString(parseFormula(f)).replace(" ","").replace("\n","").replace("/","")
                    counts = [len(re.findall(x,treeRep)) for x in feature_list]
                    go_dict_ast[el].append(counts)
                    go_dict[el].append(f)
                
            for i in range(0, model.getNumCompartments()):
                sp = model.getCompartment(i)
                if cmprt != 'all':
                    if sp.getId() in cmprt:
                        cpfor[sp.getId()].append(formulas)
                else:
                    cpfor[sp.getId()].append(formulas)

    ## get compartment|formula structure
    return (cpfor,go_dict,go_dict_ast)


def inter_component_distances(formula_file,measure="ED",precomputed=None,jid="default"):
    
    import seaborn as sns
    if precomputed == None:
        if measure == "ED":
            import editdistance as ed
        elif measure == "fuzzy" or measure == "fuzzy_plain":
            from fuzzywuzzy import fuzz
            from fuzzywuzzy import process
        else:
            print("Default measure..")

        import numpy as np

        print("Starting pairwise distance measurements..")
        distframe = pd.DataFrame()
        ## double loop for pairwise distances v is of form list of lists
        partial = 0
        totlen = len(formula_file.keys())
        print("Number of formula clusters {}".format(totlen))
        for k,v in formula_file.items():
            partial+=1
            if partial % 1 == 0:
                print(float((partial*100)/totlen),"%","complete.")
            for k2,v2 in formula_file.items():

                ## first get representative formulas for individual components
                distMinAvg = 0                            
                all_pairs = []
                for f1 in v:
                    for f2 in v2:
                        if (f1,f2) not in all_pairs:
                            all_pairs.append((f1,f2))

                if measure == "ED":
                    distMinAvg = np.mean([pool.apply(ed.eval, args=(x,y,)) for x,y in all_pairs])
                if measure == "fuzzy":
                    distMinAvg = np.mean([pool.apply(fuzz.partial_ratio, args=(x,y,)) for x,y in all_pairs])
                else:
                    distMinAvg = np.mean([pool.apply(fuzz.ratio, args=(x,y,)) for x,y in all_pairs])

                
                if distMinAvg >= 0:
                    distframe = distframe.append({'First component' : k, 'Second component' : k2, 'distance' : distMinAvg},ignore_index=True)
    else:
        distframe = pd.read_csv(precomputed)
 

    distframe.to_csv("out_files/"+measure+"_"+jid+".csv")

def draw_heatmap_basic(fname):
    from sklearn import preprocessing
    scaler = preprocessing.MinMaxScaler()
    print(fname)
    distframe = pd.read_csv(fname)
    print(distframe.describe())
    print("read file..")
    distframe[['distance']] = pd.DataFrame(scaler.fit_transform(distframe[['distance']]))
    
#    distframe['distance'] = distframe['distance']/distframe['distance'].max()
#    distframe = distframe[(distframe['distance'] >= 0)]

    indata = distframe.pivot("First component","Second component","distance")
    indata=indata.fillna(0)
    rc={'axes.labelsize': 10, 'font.size': 2, 'legend.fontsize': 100, 'axes.titlesize': 0}
    plt.rcParams.update(**rc)
    #g = sns.heatmap(dframe,cmap=px,cbar_kws={'label': 'Number of edges (log10)'},linewidths=.5)
    ax = sns.heatmap(indata,cmap="Paired",)
    ax.figure.tight_layout()
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.xlabel("Second term")
    plt.ylabel("First term")
    plt.show()

    
    ax2 = sns.clustermap(indata,cmap="Paired")    
    plt.setp(ax2.ax_heatmap.yaxis.get_majorticklabels(), rotation=30)
    plt.setp(ax2.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.xlabel("Second term")
    plt.ylabel("First term")
    plt.show()

def construct_consensus_string(termlist):

    ## for each therm, do align, do top alignment,
    ## for each consequential, do print
    ## craete consensus

    pass


def compare_consensus_sequences(lconsensus):

    ## maybe even apply the inter sequence comparison, yet no inner loop here..
    pass

def get_similarity_list(fname):
    distframe = pd.read_csv(fname)
    distframe= distframe[distframe['First component'] != distframe['Second component']].sort_values(['distance'],ascending=False)
    distframe.to_csv(fname.split(".")[0]+"similar_list.csv")
    distframe.to_latex(fname.split(".")[0]+"similar_list.tex")

def fingerprints_inter(formula_file,precomputed=None,jid="default"):

    dframe = pd.DataFrame()
    print("Comparing fingerptints")
    count = 0
    for x,y in formula_file.items():
        count+=1
        if count/len(formula_file.items()) % 0.01:
            print("Processed",count/len(formula_file.items()))
        for x2,y2 in formula_file.items():
            m1 = np.mean(y, axis=0)
            m2 = np.mean(y2, axis=0)
            dsum = np.sum(np.absolute(m2-m1))
            dframe = dframe.append({'First component' : x, 'Second component' : x2, 'distance' : dsum},ignore_index=True)

    
    dframe.to_csv("out_files/"+"_"+jid+".csv")
    
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
    parser.add_argument("--draw_hm",help="Draw a basic heatmap")
    parser.add_argument("--finger",help="sign_counts")
    parser.add_argument("--job",help="job_name") 
    args = parser.parse_args()

    if args.simstats:
        get_similarity_list(args.simstats)


    if args.draw_hm:
        draw_heatmap_basic(args.draw_hm)

    
    print("Beginning extraction..")
    
    datafolder = "data/BioModels_Database-r31_pub-sbml_files/curated"
    model_getter = model_generator(datafolder)

    ## query parts of the cell/ organism
    compartments_to_check=['cell','nucleus','plasma','nuclei','CellSurface','cytosol','vacuole','Lysosome','Mitochondria','cellsurface','Endosome']
    
    if args.stats:
        ## those are some basic numeric statistics regarding individual models
        get_basic_stats(model_getter,compartment=compartments_to_check)

    ## get both term sets..
    compartment_formulas, go_formulas, ast = getModelMath(model_getter,cmprt="all")

    import random
    subset = random.sample(compartment_formulas.keys(),10)
    print(subset)
    compartment_formulas = dict((key,value) for key, value in compartment_formulas.items() if key in subset)

    ## follow either go terms or compartment names..
    if args.goterms:
        slist = list((go_formulas.keys()))
        comp_formulas = {k : go_formulas[k] for k in slist}
        print(len(comp_formulas.keys())," Individual GO terms found for processing..")
    else:
        comp_formulas = compartment_formulas
    
    if args.interlev:

        ## this only gets the saved data, which is further      
        inter_component_distances(comp_formulas,jid=args.job)

    if args.interfuzzy:
        inter_component_distances(comp_formulas,measure="fuzzy",jid=args.job)

    if args.interfuzzybasic:
        inter_component_distances(comp_formulas,measure="fuzzy_plain",jid=args.job)

    if args.finger:
        fingerprints_inter(ast,jid=args.job)
