#!/bin/bash
#
# bash (Bourne-Again SHell) script to run the programs Gauss 03 and MDCLG coupled
# 
# ATTENTION!!! READ ABOUT THIS SCRIPT AT THE END OF THE FILE
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
# !!!!!!!!! NOTICE THAT YOU ONLY HAVE TO SET UP THE VARIABLES INSIDE THE BOX [INPUT OF THE BASH SCRIPT] !!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
#                                                                                           
#          .......................................................................................               
#          .                                                                                     .
#          .       NAME YOUR CLUSTER AND DEFINE THE POSITIONS AND VELOCITIES STATUS FLAG         .     
#          .       SET UP THE NUMBER OF ATOMS AND THE NUMBER OF TIME-STEPS FOR CURRENT RUN       .
#          .       ASK OR NOT GAUSSIAN TO DO NATURAL BOND ORBITAL ANALYSIS                       . 
#          .                                                                                     .    
#          .......................................................................................     
#
#....................................... INPUT OF THE BASH SCRIPT ...............................................
#
#   ATTENTION!!! NO SPACES ARE ALLOWED IN THE DEFINITION OF VARIABLES [e.g., cluster = NaCl is wrong]
#                                                                     [e.g., cluster=NaCl is 0k     ]
#
#   CAREFULL !!! if start=1 you must set LPOST = F AND LVEL = F IN THE FILE mdcontrl.dat
#  
#   start=0     # call mdclg
#   start=1     # start from old positions and velocities
#
#   CAREFUL!!! IF nbo=1 MAKE SURE pop=nbo (instead of pop=full) IN GAUSSIAN OPTION LINES IN MDINPUT.DAT and LNBO=.TRUE.
#  
#   nbo=0       # Do not ask for Natural Bond Orbital Analysis
#   nbo=1       # Ask for Natural Bond Orbital Analysis  
#
#   MODIFICATIONS - 2/2007: Total dipole moment of Gaussian grep
#
cluster=AzDCM16
start=1                                                                                            
natoms=111
max=1
nbo=0
gaussian=g09                     
#
#....................................... DEFINE INITIAL POSITIONS AND VELOCITIES .................................
#
echo ""
echo "Born-Oppenheimer md simulation of $cluster"
echo ""
#
if [ "$start" -eq 0 ]
then
echo "Assigning atomic positions and/or atomic velocities"
./mdclg
else
echo "Starting simulation from old positions and velocities"
fi
#
nlinesa=`expr $natoms + 2`
nlinesb=`expr $natoms + 7`
#
i=1
#
echo ""
echo "Number of steps = $max"
echo "Number of atoms = $natoms"
if [ "$nbo" -eq 1 ]
then
echo "Natural bond orbital analysis asked for in this run"
fi
echo ""
echo "Starting simulation loop"
echo ""
#
#....................................... START LOOP FOR MD SIMULATION ............................................
#
while [ "$i" -le "$max" ] ; do
   $gaussian<gaussr.com>$gaussian"out".out
   
   grep -A$nlinesa "Forces (Hartrees/Bohr)" $gaussian"out".out>gaussf.com 
   
   grep -A0  "nuclear repulsion energy" $gaussian"out".out>>gaussf.com 
   
   grep -A4  "SCF Done:" $gaussian"out".out>>gaussf.com
#Mulliken g03   
   if [ "$gaussian" == "g03" ]
   then
   grep -A$nlinesa "Mulliken atomic charges:" $gaussian"out".out>>gaussf.com 
   fi
#Mulliken g09   
   if [ "$gaussian" == "g09" ]
   then
   grep -A$nlinesa "Mulliken charges and spin densities:" $gaussian"out".out>>gaussf.com 
   fi
#NBO   
   if [ "$nbo" -eq 1 ]
   then
   grep -A$nlinesb "Summary of Natural Population Analysis:" $gaussian"out".out>>gaussf.com
   fi   
   grep -A1  "Dipole moment" $gaussian"out".out>>gaussf.com 
   ./mdclg
   echo "Molecular Dynamics step = $i out of $max"
   i=`expr $i + 1`
   wait
done
echo ""
echo "Normal termination of mdclg simulation for $cluster"
echo ""
echo "To continue this MD simulation re-define the no of time-steps in the script if needed make start=1 and re-run the script."
echo "ATTENTION!!! The input files do not have to be changed to continue the simulation correctly."
echo ""
echo "...end"
echo ""



#
#....................................... DESCRIPTION OF THE SCRIPT ...............................................
#
#   Script comments
#   ----------------
#
#                The "#" can be used to comment the script with the exception of the first line
#                For any problems with the bash scripting type [man bash] in a shell
# 
#
#   Input parameters
#   ----------------
#
#   cluster             # Define the name of your system [e.g., cluster=NaCl]
#   start=0             # Call mdclg if need to build the two initial position files and/or the velocity file
#   start=1             # Skip definition of initial positions and velocities; start from old files
#   natoms=8            # Used in the grep operations to copy the right number of lines
#   max=10              # Total number of time steps for the simulation
#   nbo=0               # Do not ask Gaussian to do Natural Bond Orbital Analysis>>>>>>>>>>>>>>>>>>>>pop=full
#   nbo=1               # Ask Gaussian to do Natural Bond Orbital Analysis>>>>>>>>>>>>>>>>>>>>>>>>>>>pop=nbo
#
#   Calculate grep lines for forces and charges
#   -------------------------------------------
#
#   nlinesa=natoms + 2   # Define number of lines to grep from number of atoms (Forces, Mulliken charges)
#   nlinesb=natoms + 7   # Define number of lines to grep from number of atoms (NBO charges)
#
#   Run Gaussian for the first time
#   -------------------------------
# 
#   g03<gaussr.com>g03out.out
#
#   Search and copy data from the output of Gaussian
#   ------------------------------------------------
#
#   For the first grep use ">"; for all the others use ">>" to append to the same file gaussf.com
#   The file gaussf.com is re-written and read by MDCLG every time-step
#   Add other grep searches if needed, for additional properties calculated by gaussian 03
#   Add the corresponding reading and writting command in the routine Gforce.f of MDCLG 
#   Alternatively simply grep additional properties to a different output file
#   In this case do not forget to use ">>" to store the properties every time-step without loosing previous data
#
#   The default "greps" are:
#  
#   grep -A$nlinesa "Forces (Hartrees/Bohr)" g03out.out>gaussf.com 
#   grep -A0  "nuclear repulsion energy" g03out.out>>gaussf.com 
#   grep -A3  "SCF Done:" g03out.out>>gaussf.com 
#   grep -A$nlinesa "Mulliken atomic charges:" g03out.out>>gaussf.com
#
#   Run the program MDCLG to advance the atomic positions and velocities
#   --------------------------------------------------------------------
#
#   ./mdclg
#
#   Repeat the process for the number of time-steps desired for the simulation
#   --------------------------------------------------------------------------
#
#   while [...] ; do
#   ...
#   done 
#
#....................................... END OF THE BASH SCRIPT ..................................................
