'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan
Date: 19/10/2022
Version = 1.3
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to download metadata in a table format for every sequence in a desired .FASTA file
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''


# Package imports
from Bio import SeqIO, Entrez
from Bio import Entrez
import time
import pandas as pd
import numpy as np


def find_qualifier_index(arr, target: str):
    """
    A function that finds the index of a specific target qualifier.
    :inputs: arr, an array containing qualifiers,
            target, a string containing the name of the target qualifier
    :returns: i, index of qualifier or None
    """

    for i in range(len(arr)):
        if len(arr[i]) > 1:
            for key, value in arr[i].items():
                if value == target:
                    return i

    return None


def get_metadata(file_directory: str, email_address: str):
    """
    A function that downloads metadata for each sequence in a .FASTA file in a .CSV tableau format.
    :inputs: file_directory, containing the directory of target .FASTA file
             email_address, your email address
    :returns: None
    """

    Entrez.api_key = "0479a2cd7f6d9e97e14c35e4d49eba018509"
    Entrez.email = email_address

    with open(file_directory) as file_handle:
        seq_record = SeqIO.parse(file_handle, "fasta")
        # Create dataframe
        df = pd.DataFrame()
                
        # Loop through records
        for seq_record in SeqIO.parse(file_handle, "fasta"):
                # Get record metadata
                
                entrez_handle = Entrez.efetch(db="nucleotide", id=seq_record.id, retmode="xml", retmax=10000)
                entrez_records = Entrez.parse(entrez_handle)
                
                # Store results as XML
                for entrez_record in entrez_records:
                    print("Retrieving metdata for " + str(seq_record.id) + "...", end="\r")
                    # Authors
                    try:
                        authors = entrez_record['GBSeq_references'][0]['GBReference_authors']
                    except:
                        authors = np.NAN

                    # Organism Name
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'], "organism")
                        if i is not None:
                            organism_name = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            organism_name = np.NAN
                    except:
                        organism_name = np.NAN

                    # Organelle
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'], "organelle")
                        if i is not None:
                            organelle = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            organelle = np.NAN
                    except:
                        organelle = np.NAN

                    # Isolation source
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'], "isolation_source")
                        if i is not None:
                            source = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            source = np.NAN
                    except:
                        source = np.NAN

                    # Taxon id
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'], "db_xref")
                        if i is not None:
                            taxon_id = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][i]['GBQualifier_value'].split(":")[-1]
                        else:
                            taxon_id = np.NAN
                    except:
                        taxon_id = v

                    # Country of origin
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'], "country")
                        if i is not None:
                            country = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            country = np.NAN
                    except:
                        country = np.NAN

                    # Lat/lon
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'], "lat_lon")
                        if i is not None:
                            lat_lon = entrez_record['GBSeq_feature-table'][0]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            lat_lon = np.NAN
                    except:
                        lat_lon = np.NAN

                    # Gene
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][2]['GBFeature_quals'], "gene")
                        if i is not None:
                            gene = entrez_record['GBSeq_feature-table'][1]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            gene = np.NAN
                    except:
                        gene = np.NAN

                    # Codon start
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][2]['GBFeature_quals'], "codon_start")
                        if i is not None:
                            codon_start = entrez_record['GBSeq_feature-table'][2]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            codon_start = np.NAN
                    except:
                        codon_start = np.NAN

                    # Translation table
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][2]['GBFeature_quals'], "transl_table")
                        if i is not None:
                            transl_table = entrez_record['GBSeq_feature-table'][2]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            transl_table = np.NAN
                    except:
                        transl_table = np.NAN

                    # Translation
                    try:
                        i = find_qualifier_index(entrez_record['GBSeq_feature-table'][2]['GBFeature_quals'], "translation")
                        if i is not None:
                            translation = entrez_record['GBSeq_feature-table'][2]['GBFeature_quals'][i]['GBQualifier_value']
                        else:
                            translation = np.NAN
                    except:
                        translation = np.NAN
                
                # Create row
                row = {'id' : seq_record.id, 'organism' : organism_name, 'authors' : authors, 'organelle' : organelle, 'isolation_source' : source, 
                        'taxon_id' : taxon_id, 'country' : country, 'lat_lon' : lat_lon, 'gene' : gene, 'codon_start' : codon_start, 'sequence' : seq_record.seq, 
                        'transl_table' : transl_table, 'translation' : translation}
                df = df.append(row, ignore_index=True)

                time.sleep(0.15)
        
        save_path_list = file_directory.split("\\")
        save_path = ""

        for i in range(len(save_path_list)):
            if i != len(save_path_list) - 1:
                save_path += save_path_list[i]
                save_path += "\\"

        df.to_csv(save_path + "!record_metadata.csv")
        print("Successfully saved metadata as: !record_metadata.csv in the: " + save_path + " directory!")