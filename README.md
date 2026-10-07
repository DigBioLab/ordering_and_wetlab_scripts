# useful Python scripts
various scripts that I have made that people may find useful for gene ordering and data analysis

You will need to have python with Bio, numpy, pandas, and dnachisel or use Stacey's .sif file (bioenv_sruge.sif)

If there are other scripts you think would be useful, ask Stacey!

Currently has:
1) single gene ordering - py
2) yeast library ordering - py
3) Sanger sequencing - ipynb

Single gene ordering:

    " * Generates an twist-ready fasta file for ordering individual constructs from a folder of PDBs and/or a concatenated FASTA file.\n"
    " * Appropriate overhangs for Golden Gate cloning into entry vector(s) of interest are added automatically.\n"
    " * Reverse translation is performed with dnachisel with constraints inherited from Ryan Kibler.\n"
    " * RECOMMENDED: check your GG assemblies at https://goldengate.neb.com/#!/\n"
    " * Wondering why the script is called John Bercow? https://www.youtube.com/watch?v=VYycQTm2HrM&ab_channel=TheSun\n"
    "\n"
    " * AVAILABLE ENTRY VECTORS:\n"
    " *** see Benchling>DBL Database>Cloning plasmids ***\n"
    " * EXAMPLE COMMAND: python DBL_order.py input.fasta -g gg_vector.fasta \n"
    " * If you are using cell free mix make sure to add --cell_free flag!"

    If you have not already, either
		a. Install a conda python 3.11 environment
			i. conda install biopython
			ii. conda install numpy
			iii. conda install pandas
			iv. pip install dnachisel
		b. Or download stacey's apptainer from sharepoint named bioenv_sruge.sif
	In your terminal run script, examples are below:
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

Sanger Sequencing:

you need your original sequences mapped to well, and a file mapping your traces to each well