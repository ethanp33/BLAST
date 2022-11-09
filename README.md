## 1. Installing Python

Firstly, to run these scripts you must install Python. You can install standalone Python. Or, if you are planning on doing lots of data science in Python, you can install Anaconda which packages Python alongside libraries useful for data science.
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
   - `examples` - This folder containins multiple files that can be used as examples to test the operation of all scripts included in this package. As a majority of these scripts can take hours to run, I thought to include the results of using these examples for convenience.
   - `blast_functions.py` - Module containing functions that run BLAST on FASTA files or FASTA files in a directory.
   - `combine_sequence_function.py` - Module containing a function to combining FASTA files in a directory.
   - `download_function.py` - Module containing a function that downloads FASTA files for species in a line seperated text file.
   - `main.py` - The **main** Python script. Run this to access all modules/scripts in a command-line format.
   - `metadata_functions.py` - Module containing functions that download metadata from the NCBI given a FASTA file containing multiple sequences.
   - `notebook.ipynb` - Python Notebook that can be opened through GitHub to preview some of the applications of this script. Mostly used this when testing and creating my scripts.
   - `README.MD` - Text file containing all of the documentation that you are reading now.
   - `remove_duplicates_function.py` - Module containing a function that deletes duplicate FASTA files in a directory given they have the same sequence and country of origin.

## 5. Running Biopython-Scripts
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
