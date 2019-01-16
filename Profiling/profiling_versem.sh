#!/bin/bash



# Outputs psstat file that contains Profiling content
python ./profiling/profiling_versem.py

# Takes in pstat file and creates 
./profiling/gprof2dot/gprof2dot.py -f pstats ./profiling/psstat | dot -Tsvg -o profiling/profiled_graph.svg



