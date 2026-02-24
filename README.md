

#PYTHIA8 STAR/sPHENIX SIMULATION WORKFLOWENVIRONMENT SETUP

#Before doing anything, you must set up the environment to ensure the correctlibraries and architecture are available (SDCC/RCF specific).


cd pythia8p3_from_milan/

setup 64bits
setup nfs4

#COMPILATION

#Compile the C++ generator code (main_detroit.cc). 

#This links against the STARoptimized Pythia libraries ($OPTSTAR).

g++ main_detroit.cc -I$OPTSTAR/include -L$OPTSTAR/lib -lpythia8 -std=c++11 -o myPythia

#RUNNING THE SIMULATION

#Execute the compiled binary. This will generate the text log files in the same directory.# Clean up previous runs if necessary

#rm output/output.log  //directory strucutre was removed

rm output.log 
rm output.root

# Run the simulation (default seed)
./myPythia

# OR run with a specific seed (e.g., 12345)
 ./myPythia 12345


#POST-PROCESSING (TEXT TO ROOT)Convert the flat ASCII log file into a structured ROOT file for analysis.

# Run the conversion macro
root -l 'log_to_root.C'

#ANALYSIS

#Load the resulting ROOT file to inspect events or run further analysis macros.

root -l output.root

# FILE DESCRIPTIONS
main_detroit.cc  : The main C++ code running PYTHIA8 generation.

log_to_root.C    : ROOT macro to convert text output to ROOT trees.

output/          : Directory containing the generated log 


files.output.root      : The final ROOT file containing the Event tree.


# SCHEDULER SUBMISSION (The Main Pipeline)

#We use submit_job.sh to manage the Star Scheduler (Condor) workflow. This runs on the scratch directory to avoid filling up home quota.

#Step A: Generate Logs (Simulation)

#Submits the Pythia binary to the grid. Generates flat ASCII .log files.

./submit_job.sh logs 100


#Input: myPythia binary.
#Output: output_*.log files in scratch.

#Step B: Convert to ROOT (Pre-processing)

#Converts the ASCII logs into ROOT Trees.

./submit_job.sh root


#Input: output_*.log
#Macro: log_to_root.C
#Output: output_*.root (Raw Event Trees).


#Step B-2: Merge ROOT files to fewer 1GB size files

./addroot_tree

#By default this job submits 30 jobs 


#Step C: Flow Analysis

#Runs the flow analysis macro on the converted ROOT files.

# Example: Submit analysis jobs
./submit_job.sh analysis [nFiles]


#Input: output_*.root
#Macro: flowanalysis.C
#Output: analysis_*.root (Histograms & TProfiles).

#POST-PROCESSING (Merging)

#Once the analysis jobs are finished (check with condor_q or job_*.out logs)
# merge the resulting histogram files into a single file.

cd /star/data05/scratch/ptribedy/pythia_pp200/
hadd -f final_analysis.root analysis_*.root


# FINAL FITTING & PLOTTING

#Calculate $vn$ using 2PC method and Template Fits from the merged file.

#Prerequisites

#Ensure the following macros are present:

#mergetemplatefit.C (Template fit logic)

#rootvnfit_pythia.C (Simple Fourier fit logic)

#print_tprofile.C (Extracts TProfile data)

#run_pythia_flow.sh (Driver script)

#Execution

chmod +x run_pythia_flow.sh
./run_pythia_flow.sh



