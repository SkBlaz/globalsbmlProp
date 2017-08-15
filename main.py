## this is the main analysis class

import sys
import os.path
from libsbml import *  ## the main library
import glob ## for file search
from collections import Counter
import pandas as pd
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

def get_basic_stats(sbmlfolder):

    ## read and print the basic stats for individual .xml entries

    candidate_models = glob.glob(sbmlfolder+"/*.xml") ## get only xml files
    all_levels = []
    df = pd.DataFrame()  ## this is the container for the basic stats
    
    for en, candidate in enumerate(candidate_models):
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
                    datapoint = {"FunctionDefinitions" : model.getNumFunctionDefinitions(),
                                 "unitDefinitions" : model.getNumUnitDefinitions(),
                                 "compartmentTypes" : model.getNumSpeciesTypes(),
                                 "speciesTypes" : model.getNumSpeciesTypes(),
                                 "compartments" : model.getNumCompartments(),
                                 "species" : model.getNumSpecies(),
                                 "parameters" : model.getNumParameters(),
                                 "initialAssignments" : model.getNumInitialAssignments(),
                                 "rules" : model.getNumRules(),
                                 "constraints" : model.getNumConstraints(),
                                 "reactions" : model.getNumReactions(),
                                 "events" : model.getNumEvents(),
                                 "Compartment" : sp.getId()}
                    df = df.append(datapoint,ignore_index=True)
                                
    print(Counter(all_levels))    
    print(df.describe()) ## this is to be further plotted
    g = sns.pairplot(df,hue="Compartment", kind="reg")
    plt.show()

if __name__ == "__main__":

    ## this main class runs the whole workflow statistics scheme

    get_basic_stats("data/BioModels_Database-r31_pub-sbml_files/curated")
