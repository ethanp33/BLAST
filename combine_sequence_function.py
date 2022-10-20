'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan
Date: 19/10/2022
Version = 1.0
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to combine multiple .FASTA files in a directory into a single a single .FASTA file
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''


# Package imports
import os
import time

def combine_seq(folder_directory: str):
    """
    A function that combines all sequence files in a directory into a single sequence file.
    :inputs: folder_directory, containing the path of the folder
    :returns: None
    """
    
    # Make folder if it doesn't exist yet
    if not os.path.exists(os.path.join(folder_directory + "\\combined\\")):
        os.makedirs(os.path.join(folder_directory + "\\combined\\"))

    combined_file = open(os.path.join(folder_directory + "\\combined\\", "combined_sequences.fna"), "w")
    for f in os.listdir(folder_directory):
        print("Currently scanning: " + str(os.path.join(folder_directory, f)), end="\r")

        if f.endswith('.fna'):
            fh = open(os.path.join(folder_directory, f))
            for line in fh:
                combined_file.write(line)
            fh.close()
    combined_file.close()
    print("Successfully combined all sequences into a file called combined_sequences.fna in the " + folder_directory + "\\combined" + " folder!")