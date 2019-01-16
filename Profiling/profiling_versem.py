# Importing necessary things to profile
import cProfile, pstats, io
from pstats import SortKey

## Connecting parent path
import os
import sys
# Entering path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Importing the main script
import versem


####  Profiling 
# ---------------------------

# Intialize Profiler
pr = cProfile.Profile()

# Start Profiling
pr.enable()


## Stuff to Profile
#------------
versem.main()
#------------


# End Profiling
pr.disable()

# ---------------------------

# Setup Output
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
ps.dump_stats('profiling/psstat')
print(s.getvalue())