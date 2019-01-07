#!/bin/bash

# Outputs psstat file that contains Profiling content
python ./Profiling/profiling_versem.py

# Takes in pstat file and creates 
./Profiling/gprof2dot/gprof2dot.py -f pstats psstat | dot -Tsvg -o Profiling/callgraph.svg



