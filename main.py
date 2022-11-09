'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan
Date: 20/10/2022
Version = 1.2
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to consolidate all scripts into a useable command line format
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''


# Package imports
import os

# Module imports
from blast_functions import blast, blast_multiple
from download_function import download
from remove_duplicates_function import remove_duplicates
from combine_sequence_function import combine_seq
from metadata_functions import find_qualifier_index, get_metadata


def main_menu():
    """
    A function that opens the main selection menu.
    :inputs: None
    :returns: None
    """

    # Print menu
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("Type the number of your desired operation and hit enter:")
    print("1 - Download: Downloads FASTA sequences of a desired gene for species in a given line seperated text file.")
    print("2 - Remove Duplicates: Removes duplicate FASTA sequences in a given directory given that the country of origin and sequence is the same. (Great to run after 1).")
    print("3 - Combine Sequences: Combines all FASTA files in a given directory into a single FASTA file. (Great to run after 2).")
    print("4 - Get Metadata: Downloads metadata in a table format for every sequence in a given .FASTA file (Great to run after 3).")
    print("5 - BLAST: Runs a or mutiple BLAST searches on a given file directory or directory.")
    print("Find more information about using these scripts at: https://github.com/ethanp33/Biopython-Scripts")

    # Input
    i = str(input("Type number here: "))
    if i == "1":
        download_menu()
    elif i == "2":
        remove_dupe_menu()
    elif i == "3":
        combine_seq_menu()
    elif i == "4":
        get_metadata_menu()
    elif i == "5":
        blast_menu()
    else:
        print("Error: You haven't typed a valid number corresponding to an operation. Please try again.")
        main_menu()

def return_to_menu():
    """
    A function that returns to the main menu.
    :inputs: None
    :returns: None
    """

    i = input("Type R to Return to the main menu. Else, type E to Exit")
    if i.upper() == "R":
        main_menu()
    elif i.upper() == "E":
        exit()
    else:
        exit()


def download_menu():
    """
    A function that opens the download menu.
    :inputs: None
    :returns: None
    """

    # Print header
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("You have selected option: 1 - Download: Downloads FASTA sequences of a desired gene for species in a given line seperated text file.")
    
    # Inputs
    download_file_dir = str(input("Please input the location/directory of your text file containing line seperated species names. \ni.e. examples\\nematode_names.txt\n"))
    if not os.path.isfile(download_file_dir):
        print("File doesn't exist or incorrect directory/location. Please try again.")
        download_menu()
    email = str(input("Please input your email address: \ni.e. epay0001@student.monash.edu\n"))
    gene = str(input("Please input the name of the gene you would like to search for. \ni.e. COI\n"))

    # Run download
    download(file_name=download_file_dir, email_address=email, gene_name=gene)
    
    # Return to menu
    return_to_menu()


def remove_dupe_menu():
    """
    A function that opens the remove duplicate sequence menu.
    :inputs: None
    :returns: None
    """
    
    # Print header
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("You have selected option: 2 - Remove Duplicates: Removes duplicate FASTA sequences in a given directory given that the country of origin and sequence is the same. (Great to run after 1).")

    # Inputs
    folder_dir = str(input("Please input the directory of the folder in which you would like duplicate FASTA sequences to be removed. \ni.e. examples\\downloads\\nematode_names\\COI\n"))
    if not os.path.exists(folder_dir):
        print("Folder doesn't exist or incorrect directory/location. Please try again.")
        remove_dupe_menu()
    email = str(input("Please input your email address: \ni.e. epay0001@student.monash.edu\n"))

    # Run remove_duplicates
    remove_duplicates(folder_directory=folder_dir, email_address=email)

    # Return to menu
    return_to_menu()

def combine_seq_menu():
    """
    A function that opens combine sequence menu.
    :inputs: None
    :returns: None
    """
    
    # Print header
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("You have selected option: 3 - Combine Sequences: Combines all FASTA files in a given directory into a single FASTA file. (Great to run after 2).")

    # Inputs
    folder_dir = str(input("Please input the directory of the folder in which you would like FASTA files to be combined. \ni.e. examples\\downloads\\nematode_names\\COI\n"))
    if not os.path.exists(folder_dir):
        print("Folder doesn't exist or incorrect directory/location. Please try again.")
        combine_seq_menu()
    
    # Run combine_seq
    combine_seq(folder_directory=folder_dir)

    # Return to menu
    return_to_menu()

def get_metadata_menu():
    """
    A function that opens the get metadata menu.
    :inputs: None
    :returns: None
    """
    
    # Print header
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("You have selected option: 4 - Get Metadata: Downloads metadata in a table format for every sequence in a given .FASTA file (Great to run after 3).")

    # Inputs
    metadata_file_dir = str(input("Please input the location/directory of your .FASTA file containing multiple sequences. \ni.e. examples\downloads\\nematode_names\COI\combined\combined_sequences.fna\n"))
    if not os.path.isfile(metadata_file_dir):
        print("File doesn't exist or incorrect directory/location. Please try again.")
        get_metadata_menu()
    email = str(input("Please input your email address: \ni.e. epay0001@student.monash.edu\n"))

    # Run get_metadata
    get_metadata(file_directory=metadata_file_dir, email_address=email)

    # Return to menu
    return_to_menu()


def blast_menu():
    """
    A function that opens the get blast menu.
    :inputs: None
    :returns: None
    """
    
    # Print header
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("You have selected option: 5 - BLAST: Runs a or mutiple BLAST searches on a given file directory or directory.")

    # Single or multiple BLAST
    i = input("Would you like to run a single BLAST or multiple BLAST? Type S for Single or M for Multiple\n")
    if i.upper() == "S":
        
        # Inputs
        blast_file_dir = str(input("Please input the location/directory of the .FASTA file you would like to BLAST. \ni.e. examples\\blast\\AB252222.1.fna\n"))
        if not os.path.isfile(blast_file_dir):
            print("File doesn't exist or incorrect directory/location. Please try again.")
            blast_menu()
        
        blast_type = str(input("Please input blast type: \ni.e. blastn, blastp, tblastn or tblastx\n"))
        if blast_type.upper() == "BLASTN":
            i = str(input("Please input whether you would like to run a megablast: \ni.e. Y for Yes or N for No\n"))
            if i.upper() == "Y":
                megablast = True
            elif i.upper() == "N":
                megablast = False
            else:
                megablast = False
        max_hits = str(input("Please input max hits: \ni.e. 5, or None for no maximum\n"))
        if max_hits.upper() == "NONE":
            max_hits = None
        else:
            max_hits = int(max_hits)
        
        # Run blast
        blast(file_name=blast_file_dir, blast_type=blast_type, max_hits=max_hits, megablast=megablast, e_value_threshold=None)
    elif i.upper() == "M":
        # Inputs
        blast_folder_dir = str(input("Please input the location/directory of the .FASTA file you would like to BLAST. \ni.e. examples\\blast\n"))
        if not os.path.exists(blast_folder_dir):
            print("Folder doesn't exist or incorrect directory/location. Please try again.")
            blast_menu()
        
        blast_type = str(input("Please input blast type: \ni.e. blastn, blastp, tblastn or tblastx\n"))
        if blast_type.upper() == "BLASTN":
            i = str(input("Please input whether you would like to run a megablast: \ni.e. Y for Yes or N for No\n"))
            if i.upper() == "Y":
                megablast = True
            elif i.upper() == "N":
                megablast = False
            else:
                megablast = False
        max_hits = str(input("Please input max hits: \ni.e. 5, or None for no maximum\n"))
        if max_hits.upper() == "NONE":
            max_hits = None
        else:
            max_hits = int(max_hits)

        # Run blast
        blast_multiple(folder_directory=blast_folder_dir, blast_type=blast_type, max_hits=max_hits, megablast=megablast, e_value_threshold=None)
    else:
        print("Invalid selection. Please try again.")
        blast_menu()
    
    # Return to menu
    return_to_menu()


def main():
    """
    The main welcome function that runs on startup.
    :inputs: None
    :returns: None
    """

    # Create directories
    new_folders = ["downloads\\", "blast\\"]
    for new_folder in new_folders:
        if not os.path.exists(new_folder):
            os.makedirs(new_folder)

    # Description
    print("\n")
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("Last updated: 20/10/2022")
    print("Version: 1.2")
    print("- By Ethan Payne | email: epay0001@student.monash.edu")
    print("-----------------------------------------------------------------------------------------------------------------------")
    print("Welcome to Ethan's collection of useful scripts!")
    print("This was created as part of my Bachelor of Applied Data Science degree (ADS3001) at Monash University during a remote internship for Manaakl Whenua Lancare Research New Zealand")

    main_menu()

    pass

# Startup
if __name__ =="__main__":
    main()
