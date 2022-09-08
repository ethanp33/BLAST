'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan
Date: 8/09/2022
Version = 1.1
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to take a sequence in FASTA format and perform a BLAST search, exporting the results in a csv format
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

# Imports
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import pandas as pd
from lxml import etree as et
import numpy as np

def blast(file_name: str, blast_type: str, max_hits: int = None, megablast: bool = False, e_value_threshold: float = None):
    """
    A function that runs a BLAST of desired type and saves the results as a csv file
    :inputs: file_name (including file extension), blast_type (one of blastn, blastp, tblastn or tblastx),
            max_hits, megablast, e_value_threshold
    :returns: None
    """ 

    # Error check
    if blast_type.lower() in ["blastn", "blastp", "blastx", "tblastn", "tblastx"]:
        
        # Open sequence file
        seq_file = next(SeqIO.parse(open(file_name), "fasta"))

        # Megablast
        if blast_type == "blastn":
            if megablast == True:
                query = NCBIWWW.qblast(blast_type, "nt", seq_file.seq, megablast=True)

        # Start query
        print("Starting BLAST query on NCBI database...")
        print("Please hold. This process can take up to 10 minutes to complete...")
        query = NCBIWWW.qblast(blast_type, "nt", seq_file.seq)
        print(query)

        # Store results as XML
        print("Storing results as an XML file...")
        with open(str(file_name).split(".")[0] + "_results.xml", "w") as save_file:
            blast_results = query.read()
            save_file.write(blast_results)

        print("Success! " + str(file_name).split(".")[0] + "_results.xml has been saved!")
        
        # Make raw dataframe
        print("Parsing XML file as dataframe...")
        df = pd.read_xml(str(file_name).split(".")[0] + "_results.xml", xpath=".//Hsp")
        
        # Edit raw dataframe
        # Rename and drop columns
        df = df.drop(["Hsp_num", "Hsp_qseq", "Hsp_query-frame", "Hsp_hit-frame"], axis=1)
        df = df.rename(columns={"Hsp_bit-score" : "bit_score", "Hsp_score" : "score", "Hsp_evalue" : "evalue", "Hsp_query-from" : "query_from", "Hsp_query-to" : "query_to",
        "Hsp_hit-from" : "hit_from", "Hsp_hit-to" : "hit_to", "Hsp_identity" : "identity", "Hsp_positive" : "positive", "Hsp_gaps" : "gaps", "Hsp_align-len" : "align_len",
        "Hsp_hseq" : "hit_seq", "Hsp_midline" : "midline"})

        # Add new columns based on Hit XML alltributes
        tree = et.parse(open(str(file_name).split(".")[0] + "_results.xml"))
        hit_id = tree.xpath(".//Hit/Hit_id/text()")
        hit_def = tree.xpath(".//Hit/Hit_def/text()")
        hit_accession = tree.xpath(".//Hit/Hit_accession/text()")
        hit_length = tree.xpath(".//Hit/Hit_len/text()")
        df["id"] = hit_id
        df["description"] = hit_def
        df["accession"] = hit_accession
        df["acc_len"] = hit_length

        # Calculating percent identity
        df["per_ident"] = np.round((df["identity"] / df["align_len"])*100, decimals=2)

        # Calculating query coverage
        df["query_cover"] = np.round((df["query_to"] - df["query_from"] + 1)/len(seq_file.seq)*100, decimals=2)

        # Change order of columns
        df = df[["id", "description", "bit_score", "score", "query_cover", "per_ident", "evalue", "query_from", "query_to", 
        "hit_from", "hit_to", "identity", "positive", "gaps", "align_len", "acc_len", "hit_seq", "midline", "accession"]]

        # Remove rows if e-value is over threshold
        if e_value_threshold is not None:
            df = df[df.evalue < e_value_threshold]
        
        # Cap rows/hits if there is a maximum
        if max_hits is not None:
            df = df.sort_values(by="evalue", ascending=False)
            df = df.head(max_hits)

        # Save editted dataframe
        print("Saving dataframe as csv file...")
        df.to_csv(str(file_name).split(".")[0] + "_dataframe.csv")

        # End
        print("Success! Dataframe has successfully been saved as " + str(file_name).split(".")[0] + "_dataframe.csv" + ". Press any key to exit function.")
    else:
        print("blast_type not one of blastn, blastp, blastx, tblastn or tblastx")
        blast(file_name=file_name, blast_type=blast_type, max_hits=max_hits, megablast=megablast, e_value_threshold=e_value_threshold)