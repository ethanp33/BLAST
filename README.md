# Getting Started
This section covers how to install Python, the necessary packages and how to install and run my Biopython-Scripts.

## 1. Installing Python
Firstly, to run these scripts you must install Python. You can install standalone Python. Or, if you are planning on doing lots of data science in Python, you can install Anaconda which packages Python alongside libraries useful for data science. For a majority of people, I would recommend installing standalone Python.
- [Standalone Python](https://www.python.org/downloads/)
- [Anaconda](https://www.anaconda.com/)

## 2. Importing Packages into Standalone Python
1. Open Command Prompt if you are a Windows user, or open Terminal if you are a Mac user.
2. You will need to install each package individually using PIP, which is an in-built package manager for Python.
3. Install the below packages by entering these commands individually:
   - `pip install biopython`
   - `pip install pandas`
   - `pip install lxml`
   - `pip install numpy`

## 3. Importing Packages into Anaconda
1. Open Anaconda Prompt.
2. You will need to install each package individually using PIP, which is an in-built package manager for Python.
3. All of the packages used in these scripts are already included in Anaconda, however it is worthwhile upgrading each package incase it is out of date or missing by entering these commands individually:
   - `pip install biopython --upgrade`
   - `pip install pandas --upgrade`
   - `pip install lxml --upgrade`
   - `pip install numpy --upgrade`
4. If any of the packages are not installed, you can install them by entering these commandy individually:
   - `pip install biopython`
   - `pip install pandas`
   - `pip install lxml`
   - `pip install numpy`

## 4. Downloading Biopython-Scripts
1. Click the green code button above, then select download ZIP. Or alternatively, click [this](https://github.com/ethanp33/Biopython-Scripts/archive/refs/heads/main.zip) to download the necessary files directly.
2. Extract the contents of the ZIP file to a location of your choice, make sure to remember where you extracted this!
3. Open the `Biopython-Scripts-main` foldder in the location in which you extracted the ZIP file.
4. Inside you will notice various files and folders:
   - `__pycache__` - Ignore this folder! There is nothing of use here.
   - `examples` This folder contains multiple files that can be used as examples to test the operation of all scripts included in this package. As a majority of these scripts can take hours to run with large data, I thought to include examples using a small number of species.
   - `results` - This folder contains the results of running these scripts on a much larger number of species. As a majority of these scripts can take hours to run, I thought to include the results of using these examples on a large dataset for convenience.
   - `blast_functions.py` - Module containing functions that run BLAST on FASTA files or FASTA files in a directory.
   - `combine_sequence_function.py` - Module containing a function to combining FASTA files in a directory.
   - `download_function.py` - Module containing a function that downloads FASTA files for species in a line seperated text file.
   - `main.py` - The **main** Python script. Run this to access all modules/scripts in a command-line format.
   - `metadata_functions.py` - Module containing functions that download metadata from the NCBI given a FASTA file containing multiple sequences.
   - `notebook.ipynb` - Python Notebook that can be opened through GitHub to preview some of the applications of this script. Mostly used this when testing and creating my scripts.
   - `README.MD` - Text file containing all of the documentation that you are reading now.
   - `remove_duplicates_function.py` - Module containing a function that deletes duplicate FASTA files in a directory given they have the same sequence and country of origin.

## 5. Running Main.py
1. If you installed Standalone Python, open Command Prompt on Windows or Terminal on Mac. If you installed Anaconda, open Anaconda Prompt.
2. Get the directory of `Biopython-Scripts-main` In file explorer, click the box and copy the directory text like so: 
   - ![image](https://user-images.githubusercontent.com/62312637/200745802-85d73f6d-bc8f-4dc2-b74c-6f8edff6bcca.png)
3. As you can see, in this case my working directory is `D:\Biopython-Scripts-main`.
4. Now in Command Prompt / Terminal / Anaconda Prompt, type: `cd "directory"` and hit enter. Where directory is your working  directory. In my case, I would type:
   - `cd "D:\Biopython-Scripts-main"`
5. Now type: `python main.py` to run the main script:
   - ![image](https://user-images.githubusercontent.com/62312637/200746815-6b1f9b64-e4dd-47a1-add5-a765f7e15e09.png)
6. You should now see a welcome message and a main menu like so:
   - ![image](https://user-images.githubusercontent.com/62312637/200746939-c5040aec-b1e6-446f-88b6-7e8f9a48d6eb.png)

# Using Biopython-Scripts
This section covers the useage and functionality of all the scripts accessible via the command line main menu by using the included example located in the `examples` folder and assumes you work your way through the below steps **in order**.

## 1. Downloading Sequences
The download function is used to download FASTA files containing sequences from a text file contaning a line seperated text file containing a list of species names. This function downloads sequences and places them inside the `downloads\species_names\gene_name` folder. For the download function, you can find the example file here: `examples\nematode_names.txt`. As you can see, it contains line-seperated names of species we would like to download the sequences of:
- ![image](https://user-images.githubusercontent.com/62312637/200751704-a4629d0f-a566-4ba6-bcfa-39969e193079.png)
1. Select option 1 in the main menu by typing `1` and hitting enter:
   - ![image](https://user-images.githubusercontent.com/62312637/200755627-114b3815-330d-4464-b230-03bd59fcc166.png)
2. Type the location of your text file containing line seperated species names. For example, use `examples\nematode_names.txt`:
   - ![image](https://user-images.githubusercontent.com/62312637/200755872-e95b29ef-6640-499f-9b0a-28a013eb7470.png)
3. Enter your email address:
   - ![image](https://user-images.githubusercontent.com/62312637/200755929-d2608155-e910-4302-8e3e-2b5b0fb81b00.png)
4. Enter the name of the gene you would like to download sequences for. For example, use `COI`:
   - ![image](https://user-images.githubusercontent.com/62312637/200756072-3cc8a291-5178-4b15-a57a-d8d52392d494.png)
5. The application will now download all sequences in FASTA file format, and stores them in the `downloads\species_names\gene_name` folder:
   - ![image](https://user-images.githubusercontent.com/62312637/200756592-4ae0e40a-231d-4e72-9a26-30d82479be66.png)

## 2. Removing Duplicate Sequences
The remove duplicate sequence function is used to remove FASTA files containing sequences in a directory given they have the same nucleotide sequence and country of origin.
1. Select option 2 in the main menu by typing `2` and hitting enter:
   - ![image](https://user-images.githubusercontent.com/62312637/200758009-f1de2108-726e-4f8d-a3e8-4982f5217063.png)
2. Type the location of the folder containing FASTA files you would like to remove the duplicates of. For example, use `downloads\nematode_names\COI`:
   - ![image](https://user-images.githubusercontent.com/62312637/200758425-7fff92c6-2d82-4c01-a694-ac75e6e83e1d.png)
3. Enter your email address:
   - ![image](https://user-images.githubusercontent.com/62312637/200758614-cbeabe41-1f11-4563-99e8-add0bfc2b23b.png)
4. The application will now remove all sequences in the `downloads\nematode_names\COI` folder given they have the same nucleotide sequence and country of origin:
   - ![image](https://user-images.githubusercontent.com/62312637/200758857-7029b0f6-134d-45c2-a774-0ad34e5aa4d0.png)

## 3. Combining Sequences
The combine sequences function is used to combine FASTA files containing sequences in a directory into a single FASTA file. The function places the combined FASTA file inside of the `downloads\species_names\gene_name\combined` folder.
1. Select option 3 in the main menu by typing `3` and hitting enter:
   - ![image](https://user-images.githubusercontent.com/62312637/200759068-0936ce2e-9b7e-4761-b3bd-e903412c6e68.png)
2. Type the location of the folder containing FASTA files you would like to combine into a single file. For example, use `downloads\nematode_names\COI`:
   - ![image](https://user-images.githubusercontent.com/62312637/200759183-18de0f7c-f5fe-4302-94e4-bce83c0105a5.png)
 3. The application will now create a folder in `downloads\species_names\gene_name\` called `combined` which contains a FASTA file containing all the sequences:
   - ![image](https://user-images.githubusercontent.com/62312637/200759373-bb7e6e22-0156-458a-b295-7b3fafac955b.png)

## 4. Retrieving Metadata
The get metadata function is used to download metadata for all sequences inside of a FASTA file in the format of a CSV file. The function places the CSV file next to the inputted FASTA file. For example, if you input `downloads\species_names\gene_name\combined\combined_sequences.fna`, it will place the CSV file in `downloads\species_names\gene_name\combined`.
1. Select option 4 in the main menu by typing `4` and hitting enter:
   - ![image](https://user-images.githubusercontent.com/62312637/200759875-cc8bd0c2-add1-44d5-a8ff-ef9b692246f7.png)
2. Type the location of the FASTA file containing multiple sequences you would like to retrieve metadata for. For example, use `downloads\nematode_names\COI\combined\combined_sequences.fna`:
   - ![image](https://user-images.githubusercontent.com/62312637/200760100-3b3a0ee2-6449-4235-869a-3942e1dbdcfa.png)
3. Enter your email address:
   - ![image](https://user-images.githubusercontent.com/62312637/200760175-155d1a08-809e-4710-afaa-e945ca92855e.png)
4. The application will now create a CSV file containing all of the metadata in the `downloads\species_names\gene_name\combined` folder:
   - ![image](https://user-images.githubusercontent.com/62312637/200760384-d368f104-a965-4787-a21f-bc97dad3301c.png)

## 5. BLAST
The BLAST function is used to run a BLAST on a single FASTA file or all FASTA files in a directory. As a warning, this function is extremely slow and I would reccommend running a single BLAST. The function places the results of the BLAST in CSV file format next to the inputted FASTA file. For example, if you input `blast\MN095883.1.fna`, it will place the CSV file in `blast`. For the BLAST functions, you can find the example files here: `examples\blast`.
1. Select option 5 in the main menu by typing `5` and hitting enter:
   - ![image](https://user-images.githubusercontent.com/62312637/200769568-d3b85f69-cd4b-456b-8d8f-a915d3bd052d.png)
2. Type either S or M to run a Single or Multiple BLAST respectively. For example, type `S`:
   - ![image](https://user-images.githubusercontent.com/62312637/200769608-53c9df4c-e4a3-4b68-ad07-9edc421b753d.png)
3. Type the location of the FASTA file you would like to BLAST. For example, type `examples\blast\MN095883.1.fna`:
   - ![image](https://user-images.githubusercontent.com/62312637/200769644-69b111fd-716d-45ed-ae86-98cfb9a6e9a3.png)
4. Enter the blast type. For example, enter: `blastn`:
   - ![image](https://user-images.githubusercontent.com/62312637/200769675-19c23993-dda0-4d65-bc9c-ce87eaa5eecf.png)
5. Enter Y or N for yes or no respectively to run a megablast. For example, enter: `N`
   - ![image](https://user-images.githubusercontent.com/62312637/200769693-a84df10a-cee5-44cf-83cd-626d2bc099bf.png)
6. Enter the maximum number of hits. For example, enter: `5`:
   - Warning: Typing None will cause this to run forever, I would avoid it if possible.
   - ![image](https://user-images.githubusercontent.com/62312637/200769717-4d9eb7bb-de9b-4412-bc16-0e67d8cd6f93.png)
7. The application will now create a CSV file containing the BLAST results in the `blast` folder:
   - ![image](https://user-images.githubusercontent.com/62312637/200770541-f39533f1-aa9b-4e8d-b8aa-b44f32f05785.png)
