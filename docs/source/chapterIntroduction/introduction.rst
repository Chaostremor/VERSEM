Overview
++++++++




The Shape Functions
-------------------

For each element, we can convert the global coordinates
:math:`\mathbf{x}` to local coordinates that range from
:math:`( \xi , \eta ) , \quad - 1 \leq \xi \leq 1,- 1 \leq \eta \leq 1`.
The conversion from local to global coordinates is governed by the
transformation using so called shape functions such that

.. math:: \mathbf { x } ( \xi , \eta ) = \sum _ { a = 1 } ^ { n _ { a } } N _ { a } ( \xi , \eta ) \mathbf { x } _ { a }.

The SEM’s signature shape functions are the Lagrange Polynomials
:math:`l^{n_l}_\alpha`:

.. math:: \ell _ { \alpha } ^ { n _ { \ell } } ( \xi ) = \frac { \left( \xi - \xi _ { 0 } \right) \cdots \left( \xi - \xi _ { \alpha - 1 } \right) \left( \xi - \xi _ { \alpha + 1 } \right) \cdots \left( \xi - \xi _ { n _ { \ell } } \right) } { \left( \xi _ { \alpha } - \xi _ { 0 } \right) \cdots \left( \xi _ { \alpha } - \xi _ { \alpha - 1 } \right) \left( \xi _ { \alpha } - \xi _ { \alpha + 1 } \right) \cdots \left( \xi _ { \alpha } - \xi _ { n _ { \ell } } \right) }

where :math:`\alpha = 0,...,n_l`, and :math:`n_l` will be 5 to 10 (in
our case).

-  using shape functions and the transformation with Jacobian to get
   from global to local integration

-  Gauss-Lobatto-Legendre Quadrature to approximate the integrals

-  

Global Matrix Assembly
----------------------

Since all elements are interconnected and dependent on each other, they
need to be assembled in a linear system of equations.

-  diagonal Mass matrix, simplifies things...

-  this, that

Then, we end up with following system of equations:

.. math:: \mathbf{M} \ddot { \mathbf{u} } + \mathbf{K} \mathbf{u} = \mathbf{f}

However, :math:`\mathbf{u}` is still dependent on time. To solve this a
simple finite difference scheme can be applied or other time marching
methods.

Time Marching
-------------

The linear system of equations can then be solved using a simple finite
difference scheme, where

.. math:: \mathbf { u } ( t + d t ) = d t ^ { 2 } \left( \mathbf { M } ^ { T } \right) ^ { - 1 } \left[ \mathbf { f } - \mathbf { K } ^ { T } \mathbf { u } \right] + 2 \mathbf { u } - \mathbf { u } ( t - d t ).

While we will focus on this scheme first of all, we want to allow for
other time schemes as well, e.g. Crank-Nicholson, to allow for more
robust time marching using other PDEs.

Parallelisation
---------------

There are two types of parallelisation that can possibly be done. The
first one is the embarassingly parallel job of integration at element
level. This should be relatively simple implement. However, the
integration at element level is not the bottleneck of the computation,
the bottleneck is the solving of the linear system when marching in
time. The second type of parallelisation, which is often used in such
cases, is domain decomposition and synchronisation. Meaning the domain
is split into subdomains which overlap and are synchronised after each
time-step.

Organization
============

Due to our decision to make this code as easy to understand and modular
as possible, we decided to write the code in Python. We will however
keep the option open to write certain ’under the hood functions’ in C or
C++, so that there is a possibility of a quick speed up.

Input
-----

The input parameter files will be in YAML format. They include
information of physics and simulation of this problem. For physics, it
contains the following information, solver type, mesh type, equation
type, Gauss - Lobatto - Legendre Points (GLL points), post-processing
software and possibly more. The mesh will be created externally.

Parameter File
^^^^^^^^^^^^^^

A YAML parameter file will be given to provide necessary settings.
Parameters will be stored as a dictionary in a configuration object,
which will include methods like get(), set(), has(), etc. Paremeters
include:

#. Compute SH or P-SV wave

#. Boundary conditions

#. Interval of output snapshots

#. Number of timesteps

#. Length of one timestep

Source File
^^^^^^^^^^^

The source file is a YAML file which contains the following information:

#. Location

#. Frequency

#. Source time function (either built-in SFT or points to an external
   STF file)

#. Starting time

#. Angle

#. Moment tensor

#. Amplification factor

Station File
^^^^^^^^^^^^

The station file is a YAML file which contains the location information.

Mesh File
^^^^^^^^^

The mesh file will be created by an external meshing software. In this
project, we will use Cubit. Cubit outputs an Exodus file (``.e``) which
will be used to create a model object

Pre-Processing
--------------

The GLL points contain the total number of points and the weight at each
point to approximate the integral with the weighted sum. We will use a
library of hard-coded GLL points and weights for the Quadrature
calculating them.

For each node, the Jacobian has to be calculated and saved since the
locations of the nodes are not changing. This is necessary to convert
the coordinates from the global coordinate system to the local
(elemental) coordinate system on each for the local GLL quadrature (as
described in section [sec:gmatassembly]).

``Model_object``
----------------

From input parameters that described the model as well as the input mesh
file are converted into a ``model_object``. A complete model consists of
an Exodus (``.e``) file, which is produced using an external mesher and
a file that describes the material. The Exodus file will be read as a
:math:`2\times N` array, (:math:`N` is the number of points, which
contains the location of each point. The material file is a
:math:`3\times N` binary array (3 parameters are P-wave velocity, S-wave
velocity and density). A model object that contains the arrays will be
created. The ``model_object`` contains methods to iterate through points
and get neighbouring points ().

Output
------

Seismograms
^^^^^^^^^^^

Seismograms record the displacement, velocity or acceleration as a function
of time. The results are stored originally in ``.npy`` files. Each file
contains a :math:`Ndim*N\times 1` array which represents the disp/vel/acc of
all dimensions at each gll point at one time (:math:`N` is the number of 
points, and :math:`Ndim` is the number of dimensions. The file name could be, 
for example, ``u0.2869035948169.npy``, where ``u`` denotes ``displacement``,
and ``0.2869035948169`` denotes ``t=0.2869035948169 s``. The corresponding
function reads the desired data type and interpolate properly to get a 
smooth seismogram against time at each user-specified location point. All 
traces are normalized and sorted by trace index. User can also specify the
start and end time for a zoomed-in inspection. Each dimension is arranged as
one column in the resulting figure. The resulting figure is stored defaultly
as a ``.pdf`` file named such as ``u_4.32-4.38s_seismo.pdf``, where
``4.32-4.38s`` denotes the start and end time. Of course user can specify
the file name of the output figure. The resulting seismogram at each point is
stored in a seperate file named like ``u_x30_y30.npy``, where ``u`` denotes
``displacement``, ``x30`` denotes coordinate in 1st dimension and ``y30`` 
denotes coordinate in 2nd dimension. This file actually contains a 
:math:`Nt\times 2` array, where :math:`Nt` is the number of points after
intepolation in time. 

Wavefield Snapshots
^^^^^^^^^^^^^^^^^^^

Snapshots of the wavefield (displacement, velocity or acceleration) can be
configured as an output as well. The cooresponding function will read the
desired data file at all gll points and do a dense spatial interpolation
in a rectangle region defined by the gll locations. Then the function
would make an animation integrating snapshots at all time steps.
Interpolation is meant to make the snapshot seem smooth. User need to 
specify the final number of points along the 1st dimension after
interpolation. User should also specify the preferred data type to read 
and prefered dimension type to show. As default, the animation would be
saved as a ``.mp4`` file named such as ``u_xdim_xnpt500.mp4``, where
``xdim`` denotes the prefered dimension type to show and ``xnpt500`` means
the final number of points along the 1st dimension after interpolation. 
User can specify the file name of the output animation as well.

Post-Processing
---------------

The post-processing software could be Paraview, Matlab or Matplotlib in
Python, which will be the default. Since the post-processing part may not
be of first priority in this project, all the resulting files or figures
are currently compatiable to Python. As a future improvement, this
post-processing section would integrate more environment.

Unit Testing
------------

The unit test of this project will be implemented with Jenkins. Everyone
will submit his own tests to the server independently, which means all
the group members will do the unit test in the sections he is
responsible for.

Schedule And Division of Labor
==============================

Schedule
--------

11 Dec, 2018: Prototype
^^^^^^^^^^^^^^^^^^^^^^^

31 Dec, 2018: Finish the main solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

9 Jan, 2019: Final version
^^^^^^^^^^^^^^^^^^^^^^^^^^

Division of Labor
-----------------

| **Congyue Cui:** Coding of the input and output design
| **Chao Song:** Coding of the visualisation and details
| **Fan Wu:** Coding of the main solver
| **Lucas Sawade:** Coding of the main solver
| **Srijan Bharati Das:** Coding of the main solver
| **ALL:** Skeleton design; unit test; Documentation

Building the Automatic Documentation with Sphinx
------------------------------------------------


