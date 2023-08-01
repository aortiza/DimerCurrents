# Dimer Currents

In this project, I'm trying to observe the difference in the current transients and steady states when there is a single defect vs a clean system. A single defect is ususally created by having an odd number of particles. 
The contents of the `Scripts` folder is the following:

* Script1: Here I check the old data to understand what was done previously. The old data contained discrete ramps of the frequency which made it unsuitable for studying the transients because the previous final state corresponded to the new initial state. 
* Script2: Test1, In this script I run a short test simulation to ensure that everything is working correctly. 
* Script3: PlainRun, This script produces the first dataset: Run1-PlainRun. It takes about an hour to run in 15 cores.
* Script4: OddRun, This script produces the second dataset: Run2-OddRun. This has an odd number of particles which ensures a defect. It also takes about an hour to run in 15 cores.
* Script5: Analyze1: This script reads the `.lammpstrj` files, calculates the instantaneous current of each particle, and stores a new trayectory insid an `.hdf` file. It also stores the index and the boundaries. It produces a large file called 
    ../Data/Processed/23_07_25_Script5_currents.hdf
It then calculates the mean values of the current for each timestep, and stores it in
    ../Data/Processed/23_07_25_Script5_current_vs_time.hdf
* Script6: Analyze2: This script takes the integrated current vs time and creates a plot. 
* Script7: Test2, In this script I run a short test simulation where several experiments are realised for each tilt. 
* Script8: PlainRun-Stat, Same as Script 2, but runs several seeds for each tilt angle. The results are stored in 
    Dynamic/DimerCurrents/Run3-PlainRun-Stat
* Script9: OddRun-Stat, Same as Script 3, but runs several seeds for each tilt angle. The results are stored in 
    Dynamic/DimerCurrents/Run4-OddRun-Stat
* Script10: Analyze3: This script is analogous to script5, but some changes must be made to account for the new variable: "exp_n". It produces two datasets:
    Dynamic/DimerCurrents/Processed/23_08_01_Script10_currents.hdf
    Dynamic/DimerCurrents/Processed/23_08_01_Script10_currents_vs_time

