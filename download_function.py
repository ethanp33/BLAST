'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan
Date: 19/10/2022
Version = 1.1
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to take a text file containing species names in line seperated format and downloads any relevant nucleotide sequences in FASTA format 
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''


# Package imports
from Bio import Entrez, SeqIO
import os
import time


def download(file_name: str, email_address: str, gene_name: str):
    """"
    A function that searches the genbank gene database and downloads multiple fasta sequence files.
    :inputs: file_name (including file extension) of a text file containing line separated species names to download sequences for,
            email_address so NCBI knows who is requesting files data if anything goes awry,
            gene_name containing the gene of interest to search for (i.e. COI)
    :returns: None
    """

    # Parameters
    Entrez.api_key = "0479a2cd7f6d9e97e14c35e4d49eba018509"
    Entrez.email = email_address
    
    # Create folder if it doesn't exist yet
    if not os.path.exists(os.path.join("downloads\\" + file_name.split(".")[0] + "\\" + gene_name + "\\")):
            os.makedirs(os.path.join("downloads\\" + (file_name.split("\\")[-1]).split(".")[0] + "\\" + gene_name + "\\"))

    # Read text file
    with open(file_name) as f:
        species_names = [line.rstrip('\n') for line in f]

    # Start search
    print("Starting search. Please hold as only 10 API requests can be made per second...")
    record_ids = []
    record_counts = 0
    for species_name in species_names:
        handle = Entrez.esearch(db="nucleotide", term=species_name + "[Orgn] AND COI[Gene]", idtype="acc", retmax=10000)
        record = Entrez.read(handle)
        record_counts += int(record["Count"])
        record_ids += record["IdList"]
        print("Found " + str(record_counts) + " records so far...", end="\r")
        handle.close()

        # Wait
        time.sleep(0.15)
    print("Search success! For a search of " + str(len(species_names)) + " species, a total of " + str(record_counts) + " genbank records containing COI genes were found!")

    # Start download
    print("Starting sequence download process. Please hold as only 10 API requests can be made per second...")
    count = 0
    for record_id in record_ids:
        if not os.path.isfile("downloads\\" + (file_name.split("\\")[-1]).split(".")[0] + "\\" + gene_name + "\\" + str(record_id) + ".fna"):
            net_handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="fasta", retmode="text", retmax=10000)
            out_handle = open("downloads\\" + (file_name.split("\\")[-1]).split(".")[0] + "\\" + gene_name + "\\" + str(record_id) + ".fna", "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            count += 1
            print("Downloading record " + str(record_id) + " Downloaded " + str(count) + " sequences so far...", end="\r")

            # Wait
            time.sleep(0.15)

    # End
    print("Success! You can find your downloaded sequences in the downloads folder where this script is located!")