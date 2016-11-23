# Chrom3D
![3D_genome_structure](http://collaslab.org/wp-content/uploads/2016/11/3D_model_illustration_hic.png)

Chrom3D is a computational framework for efficient reconstruction of 3D genome structures using chromosome contact data (Hi-C, TCC and 5C data) and optionally lamin ChIP-seq data.  

  
  
# Installation instructions
Chrom3D is tested on MacOS & Linux. If you try to install it under Windows, please send us a few words on how it went.  
  
  
## Dependencies: 
Chrom3D is dependent on the **boost library** (>=1.54). If you do not have it,  you can install it using (e.g. on Ubuntu):

`sudo apt-get install libboost-dev`
  
Otherwise, go to http://www.boost.org/users/download/
  
If you have boost installed, but still have installation problems when running `make`, find the boost directory using following command:

`whereis boost`
  
Then add the boost path to the Makefile:
  
all:
    g++ -I [insert path from 'whereis boost' here] -I $(TCLAPINCLUDE) $(SRCPATH)/Util.cpp $(SRCPATH)/Bead.cpp $(SRCPATH)/Chromosome.cpp $(SRCPATH)/Randomizer.cpp $(SRCPATH)/Constraint.cpp $(SRCPATH)/Model.cpp $(SRCPATH)/MCMC.cpp $(SRCPATH)/Chrom3D.cpp -o Chrom3D
  
Then, `make` should work.
  
  
## Software for visualization of final structures
The final structures from Chrom3D are saved in "cmm" format that can be visualized in Chimera (https://www.cgl.ucsf.edu/chimera/). R programming language (https://www.r-project.org/) can also be used to visualize the structures and preform statistical calculations, but this requires customized code (we provide examples below).  
  
  
## Installation
Download Chrom3D source package by running following command in a terminal:

`git clone https://github.com/CollasLab/Chrom3D.git`

or download the zip file from [text link] https://github.com/CollasLab/Chrom3D/archive/master.zip and extract the file to a prefered directory.
  
Enter Chrom3D directory and compile the program by typing `make` in the terminal.
  
The resulting Chrom3D binary should then be created for you in this folder.
  
Errors?  Please contact us (monikasekelja @ gmail.com).  
  
  
## Test the installation by running a test simulation
To make sure that Chrom3D is working properly on your machine, run your first simulation using following command:

`./Chrom3D -o ./test_files/testStructure.cmm ./test_files/testStructure.gtrack -r 5.0 -n 10000 -l 1000 -T 1`
  
The output should be as following (note that variations in scores can occur due to different versions of random number generator):
```
# Chr1 done!
# Chr2 done!
# beads: 200
# interactions: 138
# interactions with given distance: 0
# non-interactions with given minimum distance: 0
# lamin beads: 60
# non-lamin beads: 140
0 520.737 0 343.481 0 0 0 864.219
1000 195.47 0 7.38395 0 0 0 202.854
2000 164.545 0 6.66897 0 0 0 171.214
3000 159.279 0 5.13714 0 0 0 164.416
4000 152.387 0 4.09988 0 0 0 156.487
5000 153.899 0 6.24366 0 0 0 160.142
6000 128.948 0 5.30592 0 0 0 134.254
7000 120.901 0 6.26731 0 0 0 127.169
8000 118.694 0 5.42927 0 0 0 124.123
9000 110.4 0 5.55677 0 0 0 115.957
10000 111.08 0 7.04636 0 0 0 118.126
```  
  
   
  
# Getting started using examples
  
## Global modelling: Reconstruction of toy-genome from toy gtrack file (15 min)
* As input file, we will use toy-global-example.gtrack file found in Chrom3D/test_files directory consisting of 2x2 (diploid cell) chromosomes with information on pairwise chromsomal interactions and interactions with lamin-B-proteins
  
* To reconstruct this toy-genome, run following command (explanation below):

`./Chrom3D -o ./test_files/toy-global-example.cmm -r 3.0 --smart -n 50000 -l 5000 ./test_files/toy-global-example.gtrack`

> "-o ./test_files/toy-global-example.cmm" gives instructions to save the final output structure in cmm format in Chrom3D/test_files directory. The size of the nucleus (3.0 µm) is set using "-r 3.0". In addition to interaction and lamin constraints, we wish to use the smart constraint (--smart), so that all beads which are not associated with lamin proteins are pushed towards the center of the nucleus. We instruct the Chrom3D to run 50000 iterations (-n) and to give us information on model score every 5000 iterations. The input file is "./test_files/toy-global-example.gtrack".
  
* The output from this command should look like this (note that the variations in scores can occur due to different versions of random number generators):
```
# Chr1_A done!
# Chr1_B done!
# Chr2_A done!
# Chr2_B done!
# beads: 400
# interactions: 540
# interactions with given distance: 0
# non-interactions with given minimum distance: 0
# lamin beads: 120
# non-lamin beads: 280
0 3613.43 0 2395.8 10881.4 0 0 0 0 16890.6
5000 632.082 0 27.1099 597.226 0 0 1256.42
10000 470.293 0 27.9938 487.334 0 0 985.621
15000 413.497 0 30.0887 451.784 0 0 895.37
20000 376.615 0 31.6875 432.618 0 0 840.92
25000 368.641 0 29.8059 425.026 0 0 823.472
30000 357.779 0 29.3033 413.574 0 0 800.656
35000 345.364 0 28.5572 398.715 0 0 772.636
40000 343.694 0 28.8298 384.357 0 0 756.881
45000 330.399 0 28.4021 377.478 0 0 736.279
50000 316.824 0 27.2228 373.573 0 0 717.619
```
  
* Use Chimera to visualize the final structure:

`<chimera_install_dir>/UCSF-Chimera64-1.10.2/bin/chimera ./test_files/toy-global-example.cmm`
  
* The final structure looks like this:
[See the video of the structure](https://www.youtube.com/watch?v=VvLAZavML4Y)
  
* A tip: To visualize the transparent nuclear boundary in Chimera, go to "Favorites/Command Line" and in the command line window (bottom), write following command:

`shape sphere radius 3.0 center 0,0,0 color 1,1,0,0.4 mesh true`  
  
  
## Local modelling: Reconstruction of 3D structure of a genomic region (10 min)
* As input file, we will use GM12878.gtrack file found in Chrom3D/test_files directory. This gtrack file is constructed using 5-C data from α-globin gene locus for GM12878 cell type whose 3D conformation has previously been inferred in Bau D. et al 2011. Nat.Struct.Mol.Biiol. 
  
* To reconstruct this genomic region, run following command (explanation below).
`./Chrom3D -o ./test_files/GM12878.cmm -n 10000 -l 1000 -T 1.0 -c 0.0005 -e translate -e rotate ./test_files/GM12878.gtrack`
> The parameter "-o ./test_files/GM12878.cmm" gives instructions to save the final output structure in cmm format in the Chrom3D/test_files directory. The gtrack file contains information on interaction and non-interaction constraints, i.e. how close should a pair of interacting beads be, as well as how far should the non-interacting beads be from each other. Since we are preforming local modelling on one small chromosomal part, there are no other chromosomes or nuclear boundary to relate to so we can exclude rotate and translate move (nothing will change if we rotate or move this structure). The total number of iterations is set to 10000 iterations (-n) with informational output after each new 1000 iterations (-l). Although it is not necessary to change default temperature (-T) and cooling raten (-c) parameter values, in many cases a proper cooling rate can speed up the optimization process. To learn more about temperature and cooling rate, please go to Temperature parameter and cooling rate.
  
* The output from this command should look like this
```
# chr16 done!
# beads: 70
# interactions: 0
# interactions with given distance: 1048
# non-interactions with given minimum distance: 8076
# lamin beads: 0
# non-lamin beads: 70
0 0 0 0 0 6.05261e+08
2000 0 0 0 0 3.5495e+08
4000 0 0 0 0 3.38963e+08
6000 0 0 0 0 3.35464e+08
8000 0 0 0 0 3.3418e+08
10000 0 0 0 0 3.32864e+08
12000 0 0 0 0 3.31548e+08
14000 0 0 0 0 3.31147e+08
16000 0 0 0 0 3.29966e+08
18000 0 0 0 0 3.29339e+08
20000 0 0 0 0 3.28557e+08
```
  
* Use Chimera to visualize the final structure:

`<chimera_install_dir>/UCSF-Chimera64-1.10.2/bin/chimera ./test_files/GM12878.cmm`
  
* And the final structure visualized using [Chimera](http://www.cgl.ucsf.edu/chimera/) looks like this:
[See the video of the structure](https://www.youtube.com/watch?v=6qVumoJeqn8)  
  
  
  
## Getting started with your own data
Genome 3D reconstruction from Hi-C contact maps and lamin ChiP-seq data. Note: this section will be updated with more detailed information soon.
  
* Pre-process your Hi-C data into Hi-C contact maps
> Several great tools exist, such as https://github.com/nservant/HiC-Pro and https://github.com/MWSchmid/HiCdat
  
* Pre-process lamin ChiP-seq data
> The recommended pipeline for ChiP-seq data preprocessing and peak calling is https://github.com/CollasLab/edd
  
* Fuse Hi-C contact maps and ChiP-seq data into a gtrack file
> We will provide a python script for this soon. The gtrack file format is described in more detail [here](https://hyperbrowser.uio.no/hb/u/hb-superuser/p/gtrack/).

## Parameter information

### Constraints

* **Center constraint** is used when we wish to push all beads towards the center of the nucleus.

* **Nucleus constraint** is applied on a bead, it will seek to place itself within the nuclear boundary.

* **Smart constraint** is recommended when gtrack file contains information on lamin-associations. Then, Chrom3D will set a center constraint on all beads that are not lamin-associated.

### The random moves

The chromosome/bead moves used by Chrom3D are illustrated below:
![ChromosomeMoves](http://collaslab.org/wp-content/uploads/2016/11/RandomMoves.png)

### Temperature parameter and cooling rate
The default value of temperature T parameter is 1 and cooling rate c is 0. This means that the default optimization is based on Metropolis-Hastings algorithm since T does not influence the acceptance probability of the “bad” move. If you wish to use [simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing) approach, you need to set the value of the cooling rate c between 0 and 1. The temperature T parameter will decrease gradually during the simulation with the slope determined by c (Tnew = Tprevious*(1-c) for each accepted move. The idea behind cooling down the temperature during optimization is to speed up the convergence. It is, however, difficult to determine which cooling rate is appropriate to use and we recommend trying out different values. The rate will be dependent on the total number of iterations you wish to run during simulation and how steep you wish to narrow down the search.

We recommend to have a smooth decrease of T, obtained by setting the cooling rate to c=5/n, where n is the total number of iterations. For 100000 iterations, the cooling rate will be 0.00005 with following steepness:
![Temp1](http://collaslab.org/wp-content/uploads/2016/11/Temp1.png)

When x > 5 in c=x/n, the steepness will increase. For instance, if we choose to use c=15/n for 100000 iterations c is 0.00015, our optimization will not be able to explore the solution space very efficiently because it will be stuck in a minimum:
![Temp2](http://collaslab.org/wp-content/uploads/2016/11/Temp2.png)

When x < 5 in c=x/n, the steepness will decrease. For instance, if we choose to use c=1/n for 100000 iterations c is 0.00001, our optimization will explore the solution space in the beginning and slowly narrow down the search space:
![Temp3](http://collaslab.org/wp-content/uploads/2016/11/Temp3.png)


