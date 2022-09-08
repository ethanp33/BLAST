'''
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Created By: Ethan Payne
Date: 25/08/2022
Version = 1.0
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A Python script designed to take a sequence in FASTA format and perform a BLAST search, downloading the results in a csv format
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

# Imports
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# Open sequence file
#def open_seq():
    #try:
file_name = str(input("Enter the file name (including extension) of your sequence file. Ensure that the file is located in the same folder as this script.\n"))
seq_file = next(SeqIO.parse(open(file_name), "fasta"))
    #except 


# Blast type
def bt():
    blast_type = str(input("What type of BLAST search would you like to run? blastn, blastp, tblastn or tblastx?\n")).lower()
    if blast_type not in ["blastn", "blastp", "tblastn", "tblastx"]:
        print("Error: Response is not one of: blastn, blastp, tblastn or tblastx.")
        bt()
    else:
        return blast_type
blast_type = bt()

# Megablast
def mb():
    mb_response = str(input("Would you like to run a Megablast? Type yes or no.\n")).lower()
    if mb_response not in ["yes", "no"]:
        print("Error: Response is not one of: yes, no")
        mb()
    else:
        return mb_response

# Start query
if blast_type == "blastn":
    megablast = mb()
    if megablast == "yes":
        print("Starting query...")
        query = NCBIWWW.qblast(blast_type, "nt", seq_file.seq, megablast=True)
print("Starting query...")
query = NCBIWWW.qblast(blast_type, "nt", seq_file.seq)
print(query)

# Store results as XML
print("Storing results as an XML file...")
with open(str(file_name).split(".")[0] + "_results.xml", "w") as save_file:
    blast_results = query.read()
    save_file.write(blast_results)

print("Success! " + str(file_name).split(".")[0] + "_results.xml has been saved!")

print("Parsing XML file...")
# Parse XML results
E_VALUE_THRESH = 1e-20
for record in NCBIXML.parse(open(str(file_name).split(".")[0] + "_results.xml", "w")):
    if record.alignments:
        print("\n")
        print("query: %s" % record.query[:100])
        for align in record.alignments:
            for hsp in align.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    print("match: %s " % align.title[:100])

# Exit
i = input("Press any key to exit. ")