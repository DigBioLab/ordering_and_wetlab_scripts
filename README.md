# useful Python scripts
various scripts that I have made that people may find useful for gene ordering and data analysis

You will need to have python with Bio, numpy, pandas, and dnachisel or use Stacey's .sif file (bioenv_sruge.sif)

If there are other scripts you think would be useful, ask Stacey!

Currently has:
1) single gene ordering - py
2) yeast library ordering - py
3) Sanger sequencing - ipynb

#### Single gene ordering:

* Generates an twist-ready fasta file for ordering individual constructs from a folder of PDBs and/or a concatenated FASTA file.
* Appropriate overhangs for Golden Gate cloning into entry vector(s) of interest are added automatically.
* Reverse translation is performed with dnachisel with constraints inherited from Ryan Kibler.
* RECOMMENDED: check your GG assemblies at https://goldengate.neb.com/#!/

   1) Make a .fasta file of your amino acid sequences (if you don't have one you can easily ask claude to make you as script to convert folders of pdbs into fasta files)
	2) Download fasta of the vector you want to clone into
		a. This can be done either on Benchling
		b. Or from  Ordering and Wetlab Scripts in the cloning_vectors_DBL folder
	3) Download DBL_order.py from Ordering and Wetlab Scripts
	4) If you have not already, either
		a. Install a conda python 3.11 environment
			i. conda install biopython
			ii. conda install numpy
			iii. conda install pandas
			iv. pip install dnachisel
		b. Or download stacey's apptainer from sharepoint named bioenv_sruge.sif
	5) In your terminal run script, examples are below:
		a. Assuming you are running from the folder with the .sif file and the python file
			i. $your.fa - is the relative path to and name of the fasta with all your sequences
			ii. $plasmid.fa - is the relative path to and name of the plasmid you want to clone into
		b. With apptainer
			i. ./bioenv_sruge.sif DBL_order.py $your.fa $plasmid.fa
		c. With your own conda environment
			i. python DBL_order.py $your.fa $plasmid.fa
		d. If you are doing cell free include flag that will co-optimise for tobacco
			i. --cell_free
			ii. i.e. python DBL_order.py $your.fa $plasmid.fa --cell_free
		e. To read all the optional variables you can always use the help option
			i. python DBL_order.py -h
	6) The script will:
		a. reverse translate your sequence, codon optimizing for e.coli (option to change this)
		b. check for repeat sequences
		c. Check for internal repeats, gc content, and other things that increase complexity
		d. Pad DNA sequences to 300bp if your sequences are less than that
	7) It will give you warnings for
		a. potential alternative start sites, 
		b. lack of e280 (Trp or Ty)
		c. Not interacting with coomassie stain
		d. if there are 4 or more amino acids in a row in your base sequence
		e. If there's anything else you want it to check, let Stacey know!
	8) Script will output:
		a. Fasta file of all your dna sequences to order
		b. Csv file for your record that contains things like:
			i. well number
			ii. the actual dna fragments you ordered
			iii. the open reading frame in your plasmid
			iv. the expected translated sequences including all tags
			v. Molecular weight of expressed protein
			vi. Extinction coefficient of expressed protein
			vii. If there's more useful metrics you want, let Stacey know!


* AVAILABLE ENTRY VECTORS:
*** see [Benchling>DBL Database>Cloning plasmids](https://benchling.com/khabj/f_/BxYl1BKHiS-cloning-plasmids-p-numbers/) ***
* EXAMPLE COMMAND: python DBL_order.py input.fasta -g gg_vector.fasta
* If you are using cell free mix make sure to add --cell_free flag!

Other notes:
    * If you are doing cell free include flag that will co-optimise for tobacco
        * --cell_free
        * i.e. python DBL_order.py $your.fa $plasmid.fa --cell_free
    * To read all the optional variables you can always use the help option
        * python DBL_order.py -h

#### Sanger Sequencing:

you need your original sequences mapped to well, and a file mapping your traces to each well
