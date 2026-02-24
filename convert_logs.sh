#!/bin/bash

# Convert all output_*.log in scratch to output_*.root using log_to_root.C
# Run this on an rcas node, NOT through SUMS.

OUT_DIR="/star/data05/scratch/ptribedy/pythia_pp200"
# Directory where this script and log_to_root.C live
SRC_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "[i] Converting logs to ROOT locally"
echo "[i] OUT_DIR = ${OUT_DIR}"
echo "[i] SRC_DIR = ${SRC_DIR}"

# STAR / ROOT env
setup 64bits
setup nfs4

cd "${OUT_DIR}" || { echo "[!] Cannot cd to ${OUT_DIR}"; exit 1; }

# Sanity check: macro exists in SRC_DIR
if [ ! -f "${SRC_DIR}/log_to_root.C" ]; then
    echo "[!] log_to_root.C not found in ${SRC_DIR}"
    echo "    Make sure log_to_root.C is in the same directory as convert_logs.sh"
    exit 1
fi

# Sanity check: any logs?
shopt -s nullglob
logs=(output_*.log)
shopt -u nullglob

if [ ${#logs[@]} -eq 0 ]; then
    echo "[!] No output_*.log files found in ${OUT_DIR}"
    exit 1
fi

echo "[i] Found ${#logs[@]} log files"
echo

# ----------------------------------------------------------------------
# 1) Test conversion on first log
# ----------------------------------------------------------------------
first_log="${logs[0]}"
test_root="test_${first_log%.log}.root"
echo "[i] Test run: ${first_log} -> ${test_root}"
root -l -b -q "${SRC_DIR}/log_to_root.C(\"${first_log}\",\"${test_root}\")"

if [ ! -f "${test_root}" ]; then
    echo "[!] Test conversion failed. Check ROOT output above."
    exit 1
else
    echo "[i] Test conversion OK, removing test file"
    rm -f "${test_root}"
fi
echo

# ----------------------------------------------------------------------
# 2) Full conversion loop
# ----------------------------------------------------------------------
for log in "${logs[@]}"; do
    rootfile="${log%.log}.root"
    echo "[i] Converting ${log} -> ${rootfile}"

    if [ -f "${rootfile}" ]; then
        echo "    [skip] ${rootfile} already exists"
        continue
    fi

    root -l -b -q "${SRC_DIR}/log_to_root.C+(\"${log}\",\"${rootfile}\")"

    if [ -f "${rootfile}" ]; then
        echo "    [ok] Created ${rootfile}"
    else
        echo "    [ERR] Failed to create ${rootfile}"
    fi
done

echo
echo "[i] Done. Current files in ${OUT_DIR}:"
ls -lh output_*.log output_*.root 2>/dev/null

