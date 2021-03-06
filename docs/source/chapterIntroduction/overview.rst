Quick-Start
===========

The following will give an overview over the main functions that 
are used for seismic wave propagation using VERSEM. More PDEs are
to be implemented in the future.


Installation and Makefile
+++++++++++++++++++++++++

The easiest way to test and install all required software, test everything 
and compile the documentation is to run

``bash$ make all``

``make all`` installs all required software, unittests all tested software
in ``tests/`` and compiles the html documentation located in ``docs/`` 
(it does not compile the latexpdf automatically, our Jenkins CI server does not
provide a function that ``Sphinx`` needs, which is called ``latexmk``; hence, 
the automatic building of the documentation would fail).

There are a couple of other options using the ``Makefile``. 

- ``make init`` will install the required software using ``pip install``

- ``make docs`` will generate the automatic documentation of all modules 
  in the packages ``input``, ``src`` and ``tests`` using apidoc. The
  autogenerated ``.rst`` files are then used to compile ``html`` documentation.
  The ``html`` file is then located in the directory ``docs/build/html``.

- ``make latexpdf`` will generate a latex pdf containing the documentation.
  The document will then be located in the directory ``docs/build/latex``.

- ``make cleandocs`` cleans out both the auto-generated files as well as the
  compiled documents.

- ``make cleanresults`` cleans the results folder

- ``make clean`` cleans both results and documents.


The main driver ``versem.py``
+++++++++++++++++++++++++++++

``versem.py`` is the function that controls all computations. Each step is described 
in detail in the theory chapter.

Usage:

``bash$ python versem.py``

**Step-by-Step Walk through**

As the first step, ``versem.py`` reads the inputs located in ``src/inputs/config.json``. 

#. A new mesh is created using, an input Finite Element Mesh, an input velocity 
   model and the GLL setup. 
    
#. Following the meshing,

      #. the mass matrix (``src.mass_matrix``),
      #. stiffness matrix (``src.stiffness_matrix``) and
      #. the force vector (``src.force``)

   are constructed. All of the above take in inputs from the input file or rather the 
   in step 1 created mesh, i.e., the mass matrix needs input in form of density 
   :math:`\rho` at each GLL point on the global mesh, the stiffness matrix needs seismic 
   velocities (which are converted to Lame parameters), and the force vector needs both 
   force location as well as as the force acting at the force location. The computation
   of the stiffness matrix is the most computationally expensive part in the program 
   especially for large meshes (at the moment).

#. The nect step is actually solving the PDE in time. For this, different time schemes
   can be used as of now the default one is the Newmark time scheme, because of its 
   efficiency. More detail can be found in the time marching section in the theory 
   chapter.

When time-marching is finished, the results are saved in the subdirectory ``results/``.
In fact, three ``.npy`` files should show up and a ``timesteps/`` directory. 
``force_location.npy`` and ``fore_term.npy`` contain the force information as specified 
in the input file, ``gll_coordinates.npy`` contains a matrix with two columns specifying 
the X and Y coordinates at which the displacement measured. ``timesteps/`` contains ``.npy``
files with the displacement at each timestep. The file format is ``u<actual_simulation_time>.npy``.
Each timestep file contains a vector ``u`` with displacement of the format 
:math:`[u_{x_1},u_{y_1},u_{x_2},u_{y_2},u_{x_3},\dots]`. 

Before rerunning ``versem.py`` again, make cleanresults has to be run. Since the results will be 
saved in the same folder, previous results will not be overwritten(!). I.e., save your results to 
a new location if you want to keep them.

.. figure:: figures/VERSEM.png
      :alt: versemOverview

      A simple overview over the algorithm setup


Visualisation
+++++++++++++

There are two main formats to display the computed data. Wavefield snapshots and seismograms.
Each of them can be created using the output put from. A detailed description on how to use 
the functions can be found in the output section. See output for the plotting functions!
