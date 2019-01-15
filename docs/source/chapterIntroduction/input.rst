Input
-----

The input parameter files will be in JSON format. They include
information of physics and simulation of this problem. For physics, it
contains the following information, solver type, mesh type, equation
type, Gauss - Lobatto - Legendre Points (GLL points), post-processing
software and possibly more. The mesh will be created externally.

Configuration File
^^^^^^^^^^^^^^^^^^

A JSON parameter file will be given to provide necessary settings.
Location of the configuration file can be specified in the command line arguments. Default location is input/config.json

- ``mesh_file`` path to the mesh file

- ``material_file`` path to the material file

- ``ngll_x`` number of gll integration pointe in x direction

- ``ngll_y`` number of gll integration pointe in y direction

- ``ngll_z`` number of gll integration pointe in z direction

- ``dim`` number dimensions

- ``solver`` solver name (euler_explicit | rk2 | rk4 | ab2 | ab3 | ab4 | newmark)

- ``output_dir`` output directory

- ``boundary`` boundary type of each boundary (0: free surface | 1: internal boundary)

- ``nt`` number of timestep

- ``dt`` length of timestep (set to a negative value to compute dt automatically)

- ``sources`` defination of sources

    - ``stf_file`` python file that defines source time function

    - ``stf_args`` arguments that are passed to stf_file

    - ``location`` location of the source

    - ``term`` force direction in 2D

- ``stations`` defination of station locations

- ``trace_time_window`` time window of seismogram

- ``trace_interpolation`` interpolation factor of seismogram

- ``trace_type`` type of seismogram (u: displacement | v: velocity | a: acceleration)

- ``snap_xnpt`` number of interpolation points in x direction when plotting snapshop

- ``snap_type`` type of snapshot (u: displacement | v: velocity | a: acceleration)

- ``snap_comp`` component of snapshot (x: x direction | z: z direction | all: both)

Mesh File
^^^^^^^^^

The mesh file will be created by an external meshing software. In this
project, we will use Cubit. Cubit outputs an Exodus file (``.e``) which
will be used to create a model object.

Material File
^^^^^^^^^^^^^

The material file is a numpy array of defination of materials provided in the mesh file.

Source time function file

The source time function file is a python file which contains the function that computes the force of source over time.
Parameters that are passed to this function can be specified in the configuration file.