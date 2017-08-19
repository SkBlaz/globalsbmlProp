## this is the main analysis class

import sys
import os.path
from libsbml import *  ## the main library
import glob ## for file search
from collections import Counter
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
         print("FunctionDefinition " + str(n) + ", " + fd.getId());
 
         math = fd.getMath();
 
         # Print function arguments. 
         if (math.getNumChildren() > 1):
             print("(" + (math.getLeftChild()).getName());
 
             for n in range (1, math.getNumChildren()):
                 try:
                     print(", " + (math.getChild(n)).getName());
                 except:
                     pass
 
         print(") := ");
 
         # Print function body. 
         if (math.getNumChildren() == 0):
             print("(no body defined)");
         else:
             math = math.getChild(math.getNumChildren() - 1);
             formula = formulaToString(math);
             print(formula + "\n");

def getModelMath(genModels):

    ## this function obtaines reaction and other math related to specific models and saves them into a dataframe

    # all_levels = []
    # df = pd.DataFrame()  ## this is the container for the basic stats
    # all_compartments = []
    
    modelsave = {}
    for candidate in genModels:
        document = readSBML(candidate);
        if (document.getNumErrors() > 0):
            pass
        else:
            model = document.getModel();
            # for n in range(0, model.getNumReactions()):
            #     printReactionMath(n + 1, model.getReaction(n));
            for n in range(0,model.getNumFunctionDefinitions()):
                printFunctionDefinition(n + 1, model.getFunctionDefinition(n));
if __name__ == "__main__":


    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--stats",help="Some basic statistics")
    parser.add_argument("--math",help="Math based process clustering")    
    args = parser.parse_args()


    print(args)
    
    datafolder = "data/BioModels_Database-r31_pub-sbml_files/curated"
    model_getter = model_generator(datafolder)    
    
    if args.stats:
        ## those are some basic numeric statistics regarding individual models
        compartments_to_check=['cell','nucleus','plasma','nuclei','CellSurface','cytosol','vacuole','Lysosome','Mitochondria','cellsurface','Endosome']
        get_basic_stats(model_getter,compartment=compartments_to_check)


    ## this only gets the saved data, which is further processed.
    getModelMath(model_getter)


    ## 1.) abstract formulas
    ## 2.) formula length
    ## 3.) regex for forms and count that
    ## 4.) 
