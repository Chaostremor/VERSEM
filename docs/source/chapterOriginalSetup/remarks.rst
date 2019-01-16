Schedule And Division of Labor
==============================

Original Schedule
-----------------

11 Dec, 2018: Prototype
^^^^^^^^^^^^^^^^^^^^^^^

31 Dec, 2018: Finish the main solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We were not close to finishing the main solver

9 Jan, 2019: Final version
^^^^^^^^^^^^^^^^^^^^^^^^^^

Division of Labor
-----------------

|- **Congyue Cui:** Coding of the input and output design
|- **Chao Song:** Coding of the visualisation and details
|- **Fan Wu:** Time-scheme
|- **Lucas Sawade:** GLL-Library,  Sphinx, Jenkins & Github
|- **Srijan Bharati Das:** Construction of the Matrices
|- **ALL:** Skeleton design; unit test; Documentation


Problems
========

* A lot of the parsed inputs are not tested, meaning that if one 
  puts a string as the ngll_points there is no test of that is an 
  integer which it has to be.

* The software is still not parallelized, which was not possible 
  given the timeframe.

* For large scale computations domain decomposition is definitely
  necessary, but it would overly complicated the program and was 
  definitely not a first priority.

* We wanted to include absorbing boundary conditions. Definitely
  no time for that--, but a simple version could be created in the 
  not so far future.

* We of course, had a bunch of trouble with merging on github due 
  conflicts. That is simply unavoidable I guess, but it was frustrating
  at times.


 