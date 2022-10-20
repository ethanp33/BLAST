'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan
Date: 19/10/2022
Version = 1.2
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to delete .FASTA files in a directory given the sequence and the country of origin is the same
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

# Package imports
from Bio import SeqIO, Entrez
import time
import os


def remove_duplicates(folder_directory: str, email_address: str):
    """
    A function that scans a directory for any duplicate sequence (with same country of origin) files, deleting them.
    :inputs: folder_directory, containing the path of the folder,
             email_address, containing your email address
    :returns: None
    """
    
    Entrez.api_key = "0479a2cd7f6d9e97e14c35e4d49eba018509"
    Entrez.email = email_address

    # Loop through each file
    file_paths_to_remove = set()

    for filename_one in os.listdir(folder_directory):
        f_one = os.path.join(folder_directory, filename_one)

        # Check if the file exists
        if os.path.isfile(f_one) and filename_one.split(".")[-1] == "fna":

            # Check if file is already listed as a duplicate
            if f_one not in file_paths_to_remove:
            
                with open(f_one) as handle_one:
                    current_seq = next(SeqIO.parse(handle_one, "fasta"))
                    # Check current_seq against every other file
                    for filename_two in os.listdir(folder_directory):
                        f_two = os.path.join(folder_directory, filename_two)
                        
                        # Check if the file exists
                        if os.path.isfile(f_two) and filename_two.split(".")[-1] == "fna":

                            # Check if file is already listed as a duplicate
                            if f_two not in file_paths_to_remove:
                            
                                # Check if the file isn't the same as the one we are checking
                                if filename_one != filename_two:
                                    
                                    with open(f_two) as handle_two:
                                        check_seq = next(SeqIO.parse(handle_two, "fasta"))
                                        
                                        print("Comparing sequences between " + filename_one + " and " + filename_two, end="\r")
                                        # Check if the sequences of the fasta files are the same:                     
                                        if current_seq.seq == check_seq.seq:
                                            
                                            # Checking current_seq origin
                                            current_seq_id = current_seq.id.split(".")[0] + "." + current_seq.id.split(".")[1]
                                            entrez_handle = Entrez.efetch(db="nucleotide", id=current_seq_id, retmode="xml", retmax=10000)
                                            entrez_records = Entrez.parse(entrez_handle)
                                            for entrez_record in entrez_records:
                                                try:
                                                    current_seq_origin = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][6]['GBQualifier_value']
                                                except IndexError:
                                                    current_seq_origin = "N/A"
                                            time.sleep(0.15)

                                            # Checking check_seq origin
                                            check_seq_id = check_seq.id.split(".")[0] + "." + check_seq.id.split(".")[1]
                                            entrez_handle = Entrez.efetch(db="nucleotide", id=check_seq_id, retmode="xml", retmax=10000)
                                            entrez_records = Entrez.parse(entrez_handle)
                                            for entrez_record in entrez_records:
                                                try:
                                                    check_seq_origin = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][6]['GBQualifier_value']
                                                except IndexError:
                                                    check_seq_origin = "N/A"
                                            time.sleep(0.15)
                                            
                                            # Check if it's from the same country of origin
                                            if current_seq_origin == check_seq_origin:
                                                file_paths_to_remove.add(f_two)
            

    print("Found: " + str(len(file_paths_to_remove)) + " duplicate files to remove. Removing...")  
    
    for file_name in file_paths_to_remove:
        os.remove(file_name)
        # Wait
        time.sleep(0.15)
    print("Success!")