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

pr = cProfile.Profile()
pr.enable()

## Doing stuff
versem.main()


pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
ps.dump_stats('Profiling/psstat')
print(s.getvalue())