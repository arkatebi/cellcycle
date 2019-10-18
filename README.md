# cellcycle
Random Parametric Perturbations of Gene Regulatory Circuit Uncover State Transitions in Cell Cycle, by Ataur Katebi, Vivek Kohar, and Mingyang Lu (submitted)


RACIPE STABLE VERSION: RACIPE_3.1.1

REQUIREMENTS TO RUN RACIPE:
 1. gcc  
 2. python 3.7 or above

RACIPE INSTALLATION:
 1. git clone https://github.com/arkatebi/cellcycle
 2. Run make file: 

    make

    This will create the necessary shared object files (such as simulation_clib.so) that will be used by the RACIPE software.


CONFIGURATION FILE: racipe.cfg

   This file contains the configuration data for RACIPE. If this file is not in the working directory, the following command will create it:

   python racipe.py -M=C -I1=cellcycle.tpo

   Here, cellcycle.tpo is circuit topology from user input. In the topology file, each line specifies a regulatory interaction in the format of "A B 1". Here, A is the regulator, B is the targeted gene, and ?1? is a number specifying the type of the interaction. 1: transcriptional activation; 2: transcriptional inhibition; 3: activation by inhibiting degradation; 4: inhibition by activating degradation; 5: signaling activation; 6: signaling inhibition.

   The parameters in the configuration file can be changed according to the user need.

RACIPE COMMANDS:
1. Generate configuration file (racipe.cfg):

   python racipe.py -M=C -I1=cellcycle.tpo
2. Generate both thresholds and models (full mode):

   python racipe.py -M=A -I1=cellcycle.tpo
3. Generate and save thresholds in the .prs file:

   python racipe.py -M=T -I1=cellcycle.tpo
4. Generate models using thresholds from the .prs file

   python racipe.py -M=P -I1=cellcycle.tpo -I2=cellcycle.prs
5. Generate probabilities of the nodes (half-functional rule):

   python racipe.py -M=S -I1=cellcycle.tpo \
                         -I4=cellcycle.states.txt \
                         -I5=cellcycle.params \
                         -O=cellcycle

   Outputs are saved in cellcycle.edge.stat.txt file



OPTIONS FOR KNOCK DOWN PERTURBATIONS: 

6. Single gene knock down in (parameter range file is supplied):

   python racipe.py -M=P -I1=cellcycle.tpo -I2=cellcycle.prs -KD=Cln1

7. Double gene knock down in (parameter range file is supplied):

   python racipe.py -M=P -I1=cellcycle.tpo -I2=cellcycle.prs -KD=Cln1,Clb5

8. Single gene knock down in (full mode):

   python racipe.py -M=A -I1=cellcycle.tpo -KD=Cln1

9. Double gene knock down in (full mode):

   python racipe.py -M=A -I1=cellcycle.tpo -KD=Cln1,Clb5


OUTPUT FILES FROM RUNNING RACIPE in the full mode
1. cellcycle.states.txt - gene expression for the stable states
2. cellcycle.limitcycles.txt - gene expression for the stable limit cycles
3. cellcycle.summary.txt - summary of the number of states and limit cycles 
4. cellcycle.mpr.txt - MPR (maximum production rate) values of all 
                the genes in each model 
5. cellcycle.dnr.txt - DNR (degradation rate) values of all the genes 
                in each model
6. cellcycle.tsh.txt - TSH (threshold) values used in all the 
                interactions in each model 
7. cellcycle.hco.txt - HCO (Hill coefficient) values used in 
                all the interactions in each model 
8. cellcycle.fch.txt - FCH (fold change) values in all the 
                interactions in each model
9. cellcycle.mpr.prs, cellcycle.dnr.prs, cellcycle.tsh.prs, cellcycle.hco.prs, and cellcycle.fch.prs stores the ranges for the MPR, DNR, TSH, HCO, and FCH, respectively.


