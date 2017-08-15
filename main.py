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
#    plt.savefig('images/pairplot.png', bbox_inches='tight')

if __name__ == "__main__":

    ## this main class runs the whole workflow statistics scheme

    ## first get the model generator
    datafolder = "data/BioModels_Database-r31_pub-sbml_files/curated"
    model_getter = model_generator(datafolder) ## this is the model generator
    
    ## specify compartments of interest and plot pairs of counts!
    compartments_to_check=['cell','nucleus','plasma','nuclei','CellSurface','cytosol','vacuole','Lysosome','Mitochondria','cellsurface','Endosome']
    get_basic_stats(model_getter,compartment=compartments_to_check)
    
