# cellcycle
Modeling cellcycle gene regulatory circuit

RACIPE STABLE VERSION: RACIPE_3.1.1

REQUIREMENTS TO RUN RACIPE:
 1. gcc  
 2. python 3.7 or above

RACIPE INSTALLATION:
 1. git clone  https://github.com/arkatebi/cellcycle
 2. Run make file: 

    make

    This will create the necessary shared object files 
    (such as simulation_clib.so) that will be used by 
    the RACIPE software.

CONFIGURATION FILE: racipe.cfg
This file has the configuration information that the
RACIPE uses while running. If this file is not in the
directory, the following command will create it:
   python racipe.py -M=C -I1=cellcycle.tpo
The values in the configuration file can be changed
according to the user need.

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
   outputs are saved in cellcycle.edge.stat.txt file


KNOCK DOWN OPTIONS: 

6. Single knock down in (parameter range file is supplied):
   python racipe.py -M=P -I1=cellcycle.tpo -I2=cellcycle.prs  -KD=Cln1

7. Double knock down in (parameter range file is supplied):
   python racipe.py -M=P -I1=cellcycle.tpo -I2=cellcycle.prs  -KD=Cln1,Clb5

8. Single knock down in (full mode):
   python racipe.py -M=A -I1=cellcycle.tpo -KD=Cln1

9. Double knock down in (full mode):
   python racipe.py -M=A -I1=cellcycle.tpo -KD=Cln1,Clb5

10. Make the following changes in racipe.cfg file to run in stochastic mode: 
    (I)  ADDITTIVE_NOISE_LEVEL 30
    (II) SHOT_NOISE_LEVEL 0
    (III)INTEGRATION_METHOD:
         (1) EULER 
         (2) RUNGE_KUTTA 
             TOLERANCE 10e-12(default)
    (IV)ANNEALING_LEVEL 1 

OUTPUT FILES FROM RUNNING RACIPE S/W in full mode
1. cellcycle.states.txt - gene expression for the stable states
2. cellcycle.limitcycles.txt - gene expression for the stable limit cycles
3. cellcycle.summary.txt - summary of the number of states and limit cycles 
4. cellcycle.mpr.txt - MPR (maximum production rate) values of all 
                the genes in each models 
5. cellcycle.dnr.txt - DNR (degradation rate) values of all the genes 
                in each models
6. cellcycle.tsh.txt - TSH (threshold) values used in all the 
                interactions in each model 
7. cellcycle.hco.txt - HCO (Hill coefficient) values used in 
                all the interactions in each model 
8. cellcycle.fch.txt - FCH (fold change) values in all the 
                interactions in each model
9. cellcycle.mpr.prs, cellcycle.dnr.prs, cellcycle.tsh.prs, cellcycle.hco.prs, and 
   cellcycle.fch.prs stores the ranges for the MPR, DNR, TSH, HCO, 
   and FCH, respectively. 
