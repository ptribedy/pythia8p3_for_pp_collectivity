#!/bin/bash
# ===========================
# Backup script for pythia8p3_from_milan
# ===========================

BACKUP_DIR=/star/data01/pwg/ptribedy/pythia8p3_from_milan
DATE=$(date +%Y%m%d)
OUTFILE="backup_code_pythia8p3_from_milan_${DATE}.tar.gz"

cd "$BACKUP_DIR" || { echo "Directory not found: $BACKUP_DIR"; exit 1; }

echo "Preparing file list for backup in $BACKUP_DIR ..."

# Start fresh list
: > backup_filelist.txt

# Always include this script itself (if stored here)
SCRIPT_NAME=$(basename "$0")
[ -e "$SCRIPT_NAME" ] && echo "./${SCRIPT_NAME}" >> backup_filelist.txt

########################################
# 1) Explicit important top-level files
########################################
for f in \
  README.md \
  main_detroit.cc orig_main_detroit.cc \
  log_to_root.C orig_log_to_root.C \
  flowanalysis.C \
  12_*_flowanalysis.C \
  simple_analysis.xml work.xml submit_root.xml \
  submit_job.sh old_submit_job.sh old_working_submit_job.sh \
  old_root_part_working_submit_job.sh almost_working_submit_job.sh \
  checkwheretorunjobs.csh convert_logs.sh \
  addroot_tree old_addroot_tree addroot_tree* old_addroot_tree* \
  addme clear \
  rcas_job_load.txt file.list \
; do
  [ -e "$f" ] && echo "./$f" >> backup_filelist.txt
done

########################################
# 2) All shell / helper scripts in PWD
########################################
find . -maxdepth 1 -type f \( \
  -name "*.sh" -o -name "*.csh" -o -name "*.C" -o -name "*.cc" -o -name "*.xml" \
\) >> backup_filelist.txt

########################################
# 3) Directories with possible code
########################################
# Brian/ and junk/ : include only non-binary-like code files
for d in Brian junk; do
  if [ -d "$d" ]; then
    echo "Including non-binary files from $d/ ..."
    find "$d" -type f \
      ! -name "*.root" ! -name "*.roo" ! -name "*.log" \
      ! -name "*.zip" ! -name "*.dataset" ! -name "*.session.xml" \
      ! -name "*.so" ! -name "*.o" ! -name "*.d" \
      >> backup_filelist.txt
  fi
done

########################################
# 4) Global exclusions (safety net)
########################################
# Strip out obvious large outputs / junk
grep -Ev '\.root$|\.roo$|\.log$|\.so$|\.o$|\.d$|\.session\.xml$|\.dataset$|\.package\.zip$' \
  backup_filelist.txt > tmp && mv tmp backup_filelist.txt

########################################
# 5) Deduplicate and preview
########################################
sort -u backup_filelist.txt -o backup_filelist.txt

echo
echo "============================"
echo " Preview of files to include"
echo "============================"
head -n 20 backup_filelist.txt
echo "..."
echo "(Total files: $(wc -l < backup_filelist.txt))"
echo

########################################
# 6) Create TAR archive
########################################
echo "Creating archive: $OUTFILE ..."
tar -czf "$OUTFILE" -T backup_filelist.txt
echo "Archive created successfully."
ls -lh "$PWD/$OUTFILE"

# ANSI color codes
YELLOW='\033[1;33m'
GREEN='\033[1;32m'
RED='\033[1;31m'
CYAN='\033[1;36m'
NC='\033[0m' # No Color

echo
echo -e "${GREEN}Backup complete!${NC}"
echo -e "You can copy it from: $PWD/$OUTFILE"
echo

echo -e "${RED}To unpack safely on RCF (csh):${NC}"
echo -e "set f=\`basename $OUTFILE .tar.gz\` ; mkdir -p \\$f ; tar -xzf $OUTFILE -C \\$f"
echo
echo -e "${RED}On macOS (bash/zsh):${NC}"
echo -e 'mkdir -p "${f=$(basename '"$OUTFILE"' .tar.gz)}" && tar -xzf '"$OUTFILE"' -C "$f"'
echo
echo -e "${RED}BE CAREFUL:${NC} If you just run ${YELLOW}tar -xzf $OUTFILE${NC}"
echo -e "everything will be dumped directly into the current directory (${CYAN}$PWD${NC})."
echo

