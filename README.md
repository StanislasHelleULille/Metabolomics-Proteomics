This Python script has been used to retrieve the gene ontology for the proteins found during the anlysis by mass spectrometry. 

INSTALLATION: 
  
  unzip all the files and folders in an empty folder of your choice. Do not move copy and db folders, as the python script will fetch the input files in
  these folders. 

SYSTEM REQUIREMENT: 
  This script has only been tested on Windows 10 and Windows 11. Paths in the script being written with "\\" separators like Windows path may need to be
  changed to "/" on Linux and MacOS to work correctly. 

PYTHON REQUIREMENTS: 
  
  Python version 3.12.X or later version
  
  packages required: 
  - os
  - requests
  - time
  - biopython
  
  Please verify if these packages are installed, if not, install them using the command "pip install {package name}" (or "conda install {package name}" if you use anaconda).

CSV DATAFILE FORMAT:

  the CSV file should be in this format: 
  Protein ID1,iBAQ score 1,iBAQ score 2,[...],iBAQ score x,
  Protein ID2,iBAQ score 1,iBAQ score 2,[...],iBAQ score x,
  [...]
  Protein IDx,iBAQ score 1,iBAQ score 2,[...],iBAQ score x
  
  It can or cannot contain a header, as the header is supposedly ignored during the process, unless the header contains a protein ID present in the database. 
  
  Protein IDs should correspond to a sequence descriptions present in the database in the "db" folder, or at least to a part of the description which is unique in the database. 

DATABASE FORMAT: 

  The database should be in .fasta format

USE: 
  The database we used for retrieving the sequences in Uniprot using an automated blastP is stored in "db" folder.
  The csv file with the results obtained from proteomics analyses is stored in "sequence_list" folder
  To replicate our data, please let the database "PGSC_DM_v3.4_pep_nonredundant.fasta" inside "db" folder and the file "replicats_combine.csv" in sequence_list. 
  
  if you wish to use your own file and database, please move or delete the fasta file in the "db" folder and replace it by your own database (must be in the shape of a fasta file), and move or delete "replicats_combine.csv" from "sequence_list" folder and replace it by your own csv file (must be in coma separated csv format). 
  
  To use the script: 
  
  If you wish to only replicate the results, skip the part below:
          
          Edit the script line 12 and 13 by changing the database name and the file name you wish to use
          #Change here the database name and accession file name 
          database = directory + '/db/PGSC_DM_v3.4_pep_nonredundant.fasta'
          accessions_file = directory + '/sequence_list/replicats_combine.csv'
  
  Then, simply run it with Python using the CMD or your favourite IDE, the script with automatically begin and give you update in the CMD about the status of the process. Out of BlastP rush hour, the analysis of one containing 100 different sequences file takes between 1 and 3 minutes, and should not take more than 10 minutes. The sequence search is made 100 sequences by 100 sequences. 

