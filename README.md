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
* Wondering why the script is called John Bercow? https://www.youtube.com/watch?v=VYycQTm2HrM&ab_channel=TheSun

* AVAILABLE ENTRY VECTORS:
*** see [Benchling>DBL Database>Cloning plasmids](https://benchling.com/khabj/f_/BxYl1BKHiS-cloning-plasmids-p-numbers/) ***
* EXAMPLE COMMAND: python DBL_order.py input.fasta -g gg_vector.fasta
* If you are using cell free mix make sure to add --cell_free flag!

If you have not already, either
* Install a conda python 3.11 environment
    * conda install biopython
    * conda install numpy
    * conda install pandas
    * pip install dnachisel
* Or download stacey's apptainer from sharepoint named bioenv_sruge.sif
In your terminal run script, examples are below:
* Assuming you are running from the folder with the .sif file and the python file
    * $your.fa - is the relative path to and name of the fasta with all your sequences
    * $plasmid.fa - is the relative path to and name of the plasmid you want to clone into
    * With apptainer
        * ./bioenv_sruge.sif DBL_order.py $your.fa $plasmid.fa
    * With your own conda environment
        * python DBL_order.py $your.fa $plasmid.fa
    * If you are doing cell free include flag that will co-optimise for tobacco
        * --cell_free
        * i.e. python DBL_order.py $your.fa $plasmid.fa --cell_free
    * To read all the optional variables you can always use the help option
        * python DBL_order.py -h

#### Sanger Sequencing:

you need your original sequences mapped to well, and a file mapping your traces to each well