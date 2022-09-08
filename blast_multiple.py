'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan
Date: 8/09/2022
Version = 1.0
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to take a sequence in FASTA format and perform multiple BLAST searches, exporting the results in a csv format
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

def blast_multiple(blast_type: str, max_hits: int = None, megablast: bool = False, e_value_threshold: float = None):
    """
    A function that runs multiple BLASTs of a desired type and saves the results as a csv file.
    :inputs: folder_directory, blast_type (one of blastn, blastp, tblastn or tblastx), 
            max_hits, megablast, e_value_threshold
    :returns: None
    """ 

    # Imports
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import SeqIO
    import pandas as pd
    from lxml import etree as et
    import numpy as np
    from blast import blast

    # Get each file
    file_names = []
    count = 0
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # Check if the file exists
        if os.path.isfile(f) and filename.split(".")[-1] == "fna":
            count += 1
            print("Searching " + directory + ". Found " + str(count) + " sequences to BLAST so far.", end="\r")
            file_names.append(str(f))

    print("Please hold. Each BLAST can take up to 10 minutes to complete...")
    for file_name in file_names:
        count = 0
        print("Currenting BLASTing " + str(file_name) + ". Have BLASTed " + str(count) + " sequences so far.", end="\r")
        blast(file_name=file_name, blast_type=blast_type, max_hits=max_hits, megablast=megablast, e_value_threshold=e_value_threshold)
    print("Success! BlASTed a total of " + str(count) + " sequences.")

    return