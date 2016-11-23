# Chrom3D
Chrom3D is a computational framework for efficient reconstruction of 3D genome structures using chromosome contact data (Hi-C, TCC and 5C data) and optionally lamin ChIP-seq data.

# Installation instructions
Chrom3D is tested on MacOS & Linux. If you would like to test it under Windows, please write to us.

The final structures from Chrom3D are saved in "cmm" format that can be visualized in Chimera (https://www.cgl.ucsf.edu/chimera/). R programming language (https://www.r-project.org/) can also be used to visualize the structures and preform statistical calculations, but this requires customized code (we provide examples below).


Download Chrom3D source package by running following command in a terminal:
'git clone https://github.com/CollasLab/Chrom3D.git' 
or download the zip file from [text link] https://github.com/CollasLab/Chrom3D/archive/master.zip and extract the file to a prefered directory.

Enter Chrom3D directory and compile the program by typing 'make' in the terminal.


The resulting Chrom3D binary should then be created for you in this folder.

If you get problems with the boost library:
Chrom3D is dependent on the boost library. If you do not have it,  you can install it using (e.g. on Ubuntu):
'sudo apt-get install libboost-dev'

Otherwise, go to http://www.boost.org/users/download/

If you have boost installed, but still have installation problems with 'make':
Find the boost directory using following command:

'whereis boost'

Then add the boost path to the Makefile:

all:
    g++ -I [insert path from 'whereis boost' here] -I $(TCLAPINCLUDE) $(SRCPATH)/Util.cpp $(SRCPATH)/Bead.cpp $(SRCPATH)/Chromosome.cpp $(SRCPATH)/Randomizer.cpp $(SRCPATH)/Constraint.cpp $(SRCPATH)/Model.cpp $(SRCPATH)/MCMC.cpp $(SRCPATH)/Chrom3D.cpp -o Chrom3D

Then, 'make' should work.

# RUNNING TESTS:
Below, are two toy examples of small test-runs using Chrom3D:

Global model:
./Chrom3D -o ./test_files/toy-global-example.cmm -r 3.0 --smart -n 50000 -l 5000 ./test_files/toy-global-example.gtrack

Local model:
./Chrom3D -o ./test_files/GM12878.cmm -n 10000 -l 1000 -T 1.0 -c 0.0005 -e translate -e rotate ./test_files/GM12878.gtrack

The resulting files can be visualized using Chimera: https://www.cgl.ucsf.edu/chimera/download.html
