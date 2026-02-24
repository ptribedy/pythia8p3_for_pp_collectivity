#!/bin/bash
# P. Tribedy Dec 6,2025
# ============================================================================== 
# STAR SUBMIT SCRIPT (All-in-One)
# Usage: 
#   1. Generate Logs:    ./submit_job.sh log [N_JOBS]
#   2. Convert to Root:  ./submit_job.sh root
#   3. Run Analysis:     ./submit_job.sh analysis [FILES_PER_JOB]
# ==============================================================================

if [[ -z "$1" ]]; then
    echo "Usage: $0 [log|root|analysis] [ARG2]"

   echo "1. Generate Logs:    ./submit_job.sh log [N_JOBS] "
   echo "2. Convert to Root:  ./submit_job.sh root "
   echo "2.B You must run addroot_tree: ./addroot_tree  
   echo "3. Run Analysis:     ./submit_job.sh analysis [FILES_PER_JOB] "

    exit 1
fi

MODE="$1"
ARG2="$2"
WORK_DIR="$(pwd)"
OUT_DIR="/star/data05/scratch/ptribedy/pythia_pp200"

echo "[i] OUT_DIR: ${OUT_DIR}"

# ============================================================================== 
# MODE 1: GENERATE LOGS (RUN PYTHIA)
# ==============================================================================
if [[ "$MODE" == "log" ]]; then
    
    # Default to 10 jobs if not specified
    NJOBS=${ARG2:-10}
    
    echo "[i] Mode: LOG GENERATION"
    echo "[i] Jobs to submit: ${NJOBS}"

    # 1. Base Random Seed (Unique per submission)
    BASE_SEED=$(date +%s%N | cut -b1-9)
    echo "[i] Base Seed: ${BASE_SEED}"

    # 2. Compile Code
    echo "[i] Compiling myPythia..."
    if ! command -v root-config &> /dev/null; then
        echo "[!] Warning: ROOT environment might not be loaded (setup 64bits)."
    fi
    
    g++ main_detroit.cc -I$OPTSTAR/include -L$OPTSTAR/lib -lpythia8 -std=c++11 -o myPythia
    if [ $? -ne 0 ]; then echo "[!] Compilation failed."; exit 1; fi

    # 3. Create Directories
    mkdir -p "${OUT_DIR}/log"
    mkdir -p "${OUT_DIR}/sched"

    # 4. Generate XML
    # Using nProcesses for generation. Injecting BASE_SEED via shell variable.
    cat <<EOF > "${WORK_DIR}/simple_log.xml"
<?xml version="1.0" encoding="utf-8" ?>
<job name="PythiaGen" nProcesses="${NJOBS}" filesPerHour="5" simulateSubmission="false">
    
    <command>
        setup 64bits
        setup nfs4

        cp ${WORK_DIR}/myPythia .

        # Calculate seed: Base + JobIndex
        @ mySeed = ${BASE_SEED} + \$JOBINDEX
        echo "Running Job Index (Seed): \$mySeed"
        
        ./myPythia \$mySeed

    </command>

    <!-- Save the .log file from ./output/ to scratch -->
    <output fromScratch="*.log" toURL="file:${OUT_DIR}/" />

    <stdout URL="file:${OUT_DIR}/log/gen_${BASE_SEED}_\$JOBID.out" />
    <stderr URL="file:${OUT_DIR}/log/gen_${BASE_SEED}_\$JOBID.err" />

    <Generator>
        <Location>${OUT_DIR}/sched/</Location>
    </Generator>
</job>
EOF

    echo "[i] Submitting jobs..."
    star-submit "${WORK_DIR}/simple_log.xml"
    rm "${WORK_DIR}/simple_log.xml"
    echo "[i] Done."
    exit 0
fi

# ============================================================================== 
# MODE 2: CONVERT LOGS TO ROOT
# ==============================================================================
if [[ "$MODE" == "root" ]]; then

    echo "[i] Mode: ROOT CONVERSION"

    if [[ ! -f "${WORK_DIR}/log_to_root.C" ]]; then
        echo "[!] log_to_root.C not found!"
        exit 1
    fi

    mkdir -p "${OUT_DIR}/log"
    mkdir -p "${OUT_DIR}/sched"

    # 1. Generate File List
    LIST_FILE="${OUT_DIR}/filelist.list"
    ls -1 ${OUT_DIR}/output_*.log > ${LIST_FILE}
    NUM=$(wc -l < ${LIST_FILE})
    
    if [ "$NUM" -eq "0" ]; then
        echo "[!] No log files found in ${OUT_DIR}"
        exit 1
    fi
    echo "[i] Found ${NUM} files to process."

    # 2. Generate XML
    # Using fileListSyntax="paths" + maxFilesPerProcess="1"
    cat <<EOF > "${WORK_DIR}/simple_root.xml"
<?xml version="1.0" encoding="utf-8"?>
<job name="PythiaRoot" maxFilesPerProcess="1" fileListSyntax="paths" simulateSubmission="false">

    <command>
        setup 64bits
        setup nfs4
        
        echo "Processing: \$INPUTFILE0"

        cp ${WORK_DIR}/log_to_root.C .

        set base = \`basename \$INPUTFILE0 .log\`
        set outfile = "\${base}.root"

        echo "Output will be: \$outfile"

        # Run ROOT with escaped quotes
        root -l -b -q log_to_root.C\(\"\$INPUTFILE0\",\"\$outfile\"\)
    </command>

    <input URL="filelist:${LIST_FILE}" />

    <stdout URL="file:${OUT_DIR}/log/conv_\$JOBID.out" />
    <stderr URL="file:${OUT_DIR}/log/conv_\$JOBID.err" />

    <output fromScratch="*.root" toURL="file:${OUT_DIR}/" />

    <Generator>
        <Location>${OUT_DIR}/sched/</Location>
    </Generator>
</job>
EOF

    echo "[i] Submitting jobs..."
    star-submit "${WORK_DIR}/simple_root.xml"
    rm "${WORK_DIR}/simple_root.xml"
    echo "[i] Done."
    exit 0
fi

# ============================================================================== 
# MODE 3: ANALYSIS (Run flowanalysis.C)
# ==============================================================================
if [[ "$MODE" == "analysis" ]]; then

    echo "[i] Mode: FLOW ANALYSIS"
    
    # Files per job (default 50)
    N_PER_JOB=${ARG2:-50}
    echo "[i] Files per job: ${N_PER_JOB}"

    if [[ ! -f "${WORK_DIR}/flowanalysis.C" ]]; then
        echo "[!] flowanalysis.C not found!"
        exit 1
    fi

    mkdir -p "${OUT_DIR}/log"
    mkdir -p "${OUT_DIR}/sched"

    # 1. Generate ROOT File List (global list of all ROOT files)
    ROOT_LIST="${OUT_DIR}/rootfiles.list"
    ls -1 ${OUT_DIR}/pythia_tree_*.root > ${ROOT_LIST}
    NUM_ROOT=$(wc -l < ${ROOT_LIST})

    if [ "$NUM_ROOT" -eq "0" ]; then
        echo "[!] No ROOT files found in ${OUT_DIR}"
        exit 1
    fi
    echo "[i] Found ${NUM_ROOT} ROOT files to analyze."

    # 2. Generate XML for Analysis
    # SUMS will break ROOT_LIST into chunks of up to maxFilesPerProcess files.
    # Each job receives a local \$FILELIST containing the subset of ROOT files,
    # and flowanalysis.C already knows how to deal with .list files.
    cat <<EOF > "${WORK_DIR}/simple_analysis.xml"
<?xml version="1.0" encoding="utf-8"?>
<job name="PythiaFlow" maxFilesPerProcess="${N_PER_JOB}" fileListSyntax="paths" simulateSubmission="false">

    <command>
        setup 64bits
        setup nfs4
        
        echo "Processing File List: \$FILELIST"
        cat \$FILELIST

        cp ${WORK_DIR}/flowanalysis.C .

        # \$FILELIST is a .local.list file containing the ROOT paths for this job
        set INFILE  = "\$FILELIST"
        set OUTFILE = "analysis_\$JOBID.root"

        echo "INFILE  = \$INFILE"
        echo "OUTFILE = \$OUTFILE"

        # Run ROOT Analysis: flowanalysis(inputFiles, outputFile)
        # Use the safe quoting trick so csh doesn't choke.
        root -l -b -q 'flowanalysis.C("'"\$INFILE"'","'"\$OUTFILE"'")'
    </command>

    <!-- SUMS will read ROOT_LIST and fan it out according to maxFilesPerProcess -->
    <input URL="filelist:${ROOT_LIST}" />

    <stdout URL="file:${OUT_DIR}/log/ana_\$JOBID.out" />
    <stderr URL="file:${OUT_DIR}/log/ana_\$JOBID.err" />

    <output fromScratch="analysis_*.root" toURL="file:${OUT_DIR}/" />

    <Generator>
        <Location>${OUT_DIR}/sched/</Location>
    </Generator>
</job>
EOF

    echo "[i] Submitting jobs using simple_analysis.xml..."
    star-submit "${WORK_DIR}/simple_analysis.xml"
    # rm "${WORK_DIR}/simple_analysis.xml"
    echo "[i] Done."
    exit 0
fi

echo "[!] Unknown mode: $MODE"
echo "Usage: $0 [log|root|analysis] [ARG2]"

