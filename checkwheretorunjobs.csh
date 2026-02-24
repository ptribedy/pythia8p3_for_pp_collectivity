#!/bin/csh

# Description:
 # Checks RCAS nodes (excluding rcas6001), ranks them by job count and load average,
 # and recommends the best node to run jobs based on a weighted score.

 set tempfile = "$PWD/rcas_job_load.txt"

 # Clean up temp file
 if (-e $tempfile) then
     \rm -f $tempfile
     endif
     touch $tempfile

     echo ""
     echo "Gathering RCAS node load data..."
     echo "Node                   Jobs     LoadAvg   Score" >> $tempfile

     foreach node (`condor_status -run | grep '^rcas' | awk '{print $1}' | sort -u | grep -v rcas6001`)
         # Total jobs
             set jobs = `condor_status -run | grep $node | awk '{sum += $2} END {print sum}'`
                 if ("$jobs" == "") then
                         set jobs = 0
                             endif

                                 # Try LoadAvg, else fall back to CondorLoadAvg
                                     set load = `condor_status -long $node | grep '^LoadAvg =' | awk '{print $3}'`
                                         if ("$load" == "") then
                                                 set load = `condor_status -long $node | grep '^CondorLoadAvg =' | awk '{print $3}'`
                                                     endif
                                                         if ("$load" == "") then
                                                                 set load = 0
                                                                     endif

                                                                         # Compute score
                                                                             set score = `echo "scale=2; (0.7 * $jobs) + (0.3 * $load)" | bc`

                                                                                 printf "%-22s %-8s %-8s %-8s\n" $node $jobs $load $score >> $tempfile
                                                                                 end

                                                                                 echo ""
                                                                                 echo "Top 3 RCAS nodes to consider (excluding rcas6001):"
                                                                                 sort -k4 -n $tempfile | head -10

                                                                                 # Best node name
                                                                                 set best = `sort -k4 -n $tempfile | head -1 | awk '{print $1}'`

                                                                                 echo ""
                                                                                 echo "$best is the best choice to run your jobs."

                                                                                 # Cleanup
#                                                                                 \rm -f $tempfile

