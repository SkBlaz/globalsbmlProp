
conda install -c bioconda python-libsbml
conda install -c conda-forge fuzzywuzzy
make data
cd data; wget ftp://ftp.ebi.ac.uk/pub/databases/biomodels/releases/latest/BioModels_Database-r31_pub-sbml_files.tar.bz2
tar -xvjf BioModels_Database-r31_pub-sbml_files.tar.bz2
