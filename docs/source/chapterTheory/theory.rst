Spectral Element Method Formulation
===================================

This chapter is for documenting the theoretical aspect of the solver and
associated routines. One can refer to the equations and the
corresponding code to be clear about the implementations.

Seismic Wave Equation: FEM Solver
---------------------------------

Weak form and element matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are solving the seismic wave equation in 2D:

.. math::

   \begin{aligned}
   \rho \ddot{u}_{\rm{i}} = \partial_{\rm{j}} \tau_{\rm{ij}} + f_{\rm{i}}
   \end{aligned}

The above equation is in component form where the subscript i denotes
the component of :math:`\vec{u}`. :math:`\tau_{\rm{ij}}` is the second
order stress tensor and :math:`\vec{f}` denotes the external force
distribution. :math:`\tau_{\rm{ij}}` is expressed as:

.. math:: \tau_{\rm{ij}} = \lambda \delta_{\rm{ij}} e_{\rm{kk}} + 2 \mu e_{\rm{ij}}

and,

.. math:: e_{\rm{ij}} = \frac{1}{2}(\partial_{\rm{i}}u_{\rm{j}} + \partial_{\rm{j}}u_{\rm{i}})

Therefore, the wave equation reads as:

.. math::

   \begin{aligned}
   \rho \ddot{u}_{\rm{i}} = \partial_{\rm{j}}\bigg(\frac{\lambda \delta_{\rm{ij}}}{2}(\partial_{\rm{k}}u_{\rm{k}} + \partial_{\rm{k}}u_{\rm{k}})\bigg) + \partial_{\rm{j}}\bigg(\mu (\partial_{\rm{i}}u_{\rm{j}} + \partial_{\rm{j}}u_{\rm{i}}) \bigg) + f_{\rm{i}}\end{aligned}

For starting off, we consider constant :math:`\lambda` and :math:`\mu`.
So, we assume homogeneous and isotropic medium properties.

.. math::

   \begin{aligned}
   \rho \ddot{u}_{\rm{i}} = \lambda \partial_{\rm{i}}(\partial_{\rm{k}} u_{\rm{k}}) + \mu(\partial_{\rm{i}}(\partial_{\rm{j}}u_{\rm{j}}) + \partial_{\rm{j}}^2 u_{\rm{i}})\end{aligned}

Hereafter, we switch to the weak formulation and multiply
:math:`N_{\rm{l}}`, the :math:`l`-th shape function throughout the
equation and integrate over the area element :math:`dxdy`. We expand the
displacement :math:`\vec{u}` in terms of the components
:math:`u_{\rm{i}}(x,y) = \sum_{m=1}^{nip} u_{\rm{i}}^{\rm{m}}(x,y) N_{\rm{m}}(x,y)`,
where :math:`nip` is the total number of GLL points per element where
integration has to be carried out. Therefore, the term on the left hand
side becomes:

.. math::

   \begin{aligned}
   &=& \int N_l(x,y) \rho(x,y) \ddot{u}_{\rm{i}}(x,y)dxdy \\
   &=& \int_{-1}^{1} N_l(\xi,\eta) \rho(\xi,\eta) \sum_{m=1}^{nip} \ddot{u}_{\rm{i}}^{\rm{m}} N_{\rm{m}}(\xi,\eta) \mathcal{J}(\xi,\eta) d\xi d\eta \\
   &=& \sum_{k=1}^{nip} N_l(\xi_{\rm{k}},\eta_{\rm{k}}) \rho(\xi_{\rm{k}},\eta_{\rm{k}}) \sum_{m=1}^{nip} \ddot{u}_{\rm{i}}^{\rm{m}} N_{\rm{m}}(\xi_{\rm{k}},\eta_{\rm{k}}) \mathcal{J}(\xi_{\rm{k}},\eta_{\rm{k}}) W_{\rm{k}},\end{aligned}

where :math:`W_{\rm{k}} = w_{\rm{k1}} w_{\rm{k2}}`, the product of the
weights of the constitutive lagrange polynomials at those GLL points.
:math:`\mathcal{J}` is the determinant of the Jacobian for transforming
the integral variables from global to local coordinates. So, we sum over
all the GLL points in an element with adequate weights to perform the 2D
integral. However, the shape functions
:math:`N_{\rm{m}}(\xi_{\rm{k}},\eta_{\rm{k}}) = \delta_{\rm{mk}}` and
hence the remaining summation can be expressed as a product of a
diagonal mass matrix and a column vector :math:`u_{\rm{i}}`.

.. math:: \sum_{k=1}^{nip} \rho(\xi_{\rm{k}},\eta_{\rm{k}}) \mathcal{J}(\xi_{\rm{k}},\eta_{\rm{k}}) W_{\rm{k}} \ddot{u}_{\rm{i}}^{\rm{k}}

Hence we are finally left with the local mass matrix:
:math:`\rho(\xi_{\rm{k}},\eta_{\rm{k}}) \mathcal{J}(\xi_{\rm{k}},\eta_{\rm{k}}) W_{\rm{k}}`,
where :math:`k` is the GLL point number for an element. These will form
the diagonal elements of out matrix for each element.

The vectorial representation of the second term in 3D looks like the
following:

.. math::

   \begin{aligned}
   \int_V N_l \nabla \cdot \mathbf{\tau} dV\end{aligned}

Using Gauss’ theorem and assuming a stress free boundary condition such
that :math:`\mathbf{\tau}\cdot \mathbf{\hat{\mathbf{n}}} = 0,` where
:math:`\mathbf{\hat{\mathbf{n}}}` is the normal to the boundary surface.
So, we are left with:

.. math::

   \begin{aligned}
   \int_V N_l \nabla \cdot \mathbf{\tau} dV = -\int_V \nabla N_l \cdot \mathbf{\tau}dV\end{aligned}

In case of 2D, writing in a component form yields the equations:

.. math::

   \begin{aligned}
   &=& -\int \partial_j N_l \tau_{ij} dxdy \\
   &=& -\int \partial_j N_l(x,y) \lambda(x,y) \delta_{ij} \partial_r u_r dxdy \nonumber \\ &-& \int \partial_j N_l(x,y) \mu(x,y)(\partial_i u_j + \partial_j u_i)dxdy\end{aligned}

Let the first term be A:

.. math::

   \begin{aligned}
   &=& -\int \partial_j N_l(x,y) \lambda(x,y) \delta_{ij} \partial_r \sum_m N_m(x,y) u_r^m dxdy\end{aligned}

While the second term can be split into B and C respectively:

.. math::

   \begin{aligned}
   &-& \int \partial_j N_l(x,y) \mu(x,y) \partial_i \sum_m N_m(x,y) u_j^m \nonumber \\ &-&  \int \partial_j N_l \mu(x,y) \partial_j \sum_m N_m(x,y) u_j^m dxdy\end{aligned}

Now, we can convert the global coordinates (x,y) to local coordinate
(:math:`\xi,\eta`) using the Jacobian :math:`\mathcal{J(\xi,\eta)}` and
use the quadrature to further simplifiy the integral into:

.. math::

   \begin{aligned}
   A &= -\int_{-1}^{1} \partial_i N_l(\xi,\eta) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi,\eta} \lambda(\xi,\eta) \sum_m \partial_r N_m(\xi,\eta) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi,\eta} u_r^m \left| \mathcal{J}(\xi,\eta) \right|  d\eta d\xi \\
   &= -\sum_{k=1}^{nip} \partial_i N_l(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} \lambda(\xi_k,\eta_k) \sum_m \partial_r N_m(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} u_r^m \left| \mathcal{J} (\xi_k,\eta_k)\right| W_k \nonumber \\
   &= -\sum_{m=1}^{nip} u_{\mathbf{r}}^m \sum_{k=1}^{nip} \lambda(\xi_k,\eta_k) \partial_i N_l(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k}  \partial_{\mathbf{r}} N_m(\xi_k,\eta_k) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k}  \left| \mathcal{J}(\xi_k,\eta_k) \right| W_k \nonumber\\
   \end{aligned}

The bars around the :math:`\mathcal{J}` denote the determinant. The 
bold subscript is to remind ourselves that the index needs to be
summed over for all the components (:math:`x` and :math:`y` in 2D and
:math:`x,y` and :math:`z` in 3D). Similarly,

.. math::

   \begin{aligned}
   B &= - \int_{-1}^{1} \partial_j N_l(\xi,\eta) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi,\eta} \mu(\xi,\eta) \sum_m \partial_i N_m(\xi,\eta)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi,\eta} u_j^m \left| \mathcal{J}(\xi_k,\eta_k) \right| d\xi d\eta \\
   &= -\sum_{k=1}^{nip} \partial_j N_l(\xi_k,\eta_k) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} \mu(\xi_k,\eta_k) \sum_m \partial_i N_m(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} u_j^m \left| \mathcal{J}(\xi_k,\eta_k) \right| W_k \nonumber \\ \\
   &= -\sum_{m=1}^{nip} u_{\mathbf{j}}^m \sum_{k=1}^{nip} \mu(\xi_k,\eta_k) \partial_{\mathbf{j}} N_l(\xi_k,\eta_k) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} \partial_i N_m(\xi_k,\eta_k) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k}  \left| \mathcal{J}(\xi_k,\eta_k) \right| W_k \nonumber \\
   \end{aligned}

Finally,

.. math::

   \begin{aligned}
   C &= -\int_{-1}^{1} \partial_j N_l(\xi,\eta) \left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi,\eta} \mu(\xi,\eta) \sum_m \partial_j N_m(\xi,\eta)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi,\eta} u_i^m \left| \mathcal{J}(\xi_k,\eta_k) \right|(\xi,\eta) d\xi d\eta \\
   &= -\sum_{k=1}^{nip} \partial_j N_l(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} \mu(\xi_k,\eta_k) \sum_m \partial_j N_m(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} u_i^m \left| \mathcal{J}(\xi_k,\eta_k) \right|(\xi_k,\eta_k) W_k \nonumber\\ \\
   &= -\sum_{m=1}^{nip} u_i^m \sum_{k=1}^{nip} \mu(\xi_k,\eta_k) \partial_{\mathbf{j}} N_l(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} \partial_{\mathbf{j}} N_m(\xi_k,\eta_k)\left. \frac{\partial \mathbf{\xi}}{\partial \mathbf{x}} \right|_{\xi_k,\eta_k} \left| \mathcal{J}(\xi_k,\eta_k) \right|(\xi_k,\eta_k) W_k \nonumber\\ 
   \end{aligned}


If we notice, the expression of A, B and C have free indices as
:math:`l` and :math:`i`. Here, :math:`i` denotes the component of
:math:`\vec{u}` that we are solving for and :math:`l` denotes the row
number of the column vector :math:`u_i`.

As of know this algorith is followed strictly when computing the element
matrices.


Finding the total number of nodes to be computed:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following formula is to find the total number of computational nodes
(including the GLL points) for a rectangular mesh which have to be
solved for at each time step:

.. math::

   \begin{aligned}
   N_{total} &=& ((n_y \times N_y) - (N_y - 1)) \times (N_x + 1) \nonumber \\ &+& ((n_x - 2) \times N_x) \times (N_y + 1) \nonumber \\ + ((n_x &\times& n_y - 2 (n_y + n_x) + 4))\times (N_x \times N_y)\end{aligned}

where, :math:`n_i` is the number of nodes per element in the along i-th
coordinate and :math:`N_i` is the number of elements along the i-th
coordinate. Therefore, if we have a mesh with :math:`n_x` = :math:`n_y`
= 5 and :math:`N_x` = 20 and :math:`N_y` = 10, then :math:`N_{total}` =
3321.

Element and Node Numbering
~~~~~~~~~~~~~~~~~~~~~~~~~~

The method developed here relies on a FEM mesh created by an external
meshing software, which creates a ``.e`` Exodus file. The Exodus file
provides all the necessary information for a mesh: (a) the coordinates
of the nodes defining the control points, (b) connectivity of the
elements, as well as (c) the boundary nodes. This information is used to
create an arbitrary number of interpolation points on each element with
a new set of connectivity and numbers. This ensures that increasing
computational power can be harnessed for higher accuracy. The input mesh
therefore contains only quadrilaterals defined by four nodes (as of now;
will be extended to higher number of nodes in the future). The functions
used to read and redefine the mesh are located within the
``src/mesh_spec.py`` script. The numbering itself is performed using the
following algorithm shown in figure [fig:element\_numbering].

.. figure:: figures/element_fig_edit.pdf
   :alt: Control points of initial mesh vs GLL collocation points

   Example of GLL node and element numbering algorithm using 5x5 GLL
   points and two elements. The numbers in the red, rounded boxes denote
   the element number. The circles and white boxed numbers show the
   control point locations and numbers, respectively. the small crosses
   and numbers show the GLL point locations and numbers, respectively.

In the figure, we show the GLL points in each direction, and also the
GLL points lie on the edges as well as the interior of each element. The
locations of the GLL points contribute to the simplification of the mass
matrix, which becomes diagonal, because of the GLL quadrature used.

The interpolation of the GLL points on the Mesh of quadrilaterals is 
governed by the following equation

.. math:: \mathbf { x } ( \xi , \eta ) = \sum _ { a = 1 } ^ { n _ { a } } N _ { a } ( \xi , \eta ) \mathbf { x } _ { a }.

where :math:`N` are the basis functions, also called the shape function, 
namely, the Lagrange Polynomials. This step is applied within the module and 
function ``src.mesh_spec.mesh_interp2D()``.

Lagrange Polynomials
~~~~~~~~~~~~~~~~~~~~

The chosen test functions for the spectral element method are the
Lagrange polynomials. The Lagrange polynomials are interpolation
polynomials between certain points and defined in the following manner.

.. math::

   \label{eq:lagrange}
       \ell^n_i(\xi) = \prod_{\substack{j = 0\\ j \neq i}}^{n+1} \frac{\xi - \xi_j}{\xi_i - \xi_j},

where :math:`\xi_i` are the collocation points, :math:`i` is the number
of the polynomial with respect to collocation points and :math:`n` is
the degree of the polynomial :cite:`Komatitsch1999`. 
I.e., the elaborate expression would look like this:
The SEM’s signature shape functions are the Lagrange Polynomials
:math:`l^{n_l}_\alpha`:

.. math:: 

       \ell _ { \alpha } ^ { n _ { \ell } } ( \xi ) = \frac { \left( \xi - \xi _ { 0 } \right) \cdots \left( \xi - \xi _ { \alpha - 1 } \right) \left( \xi - \xi _ { \alpha + 1 } \right) \cdots \left( \xi - \xi _ { n _ { \ell } } \right) } { \left( \xi _ { \alpha } - \xi _ { 0 } \right) \cdots \left( \xi _ { \alpha } - \xi _ { \alpha - 1 } \right) \left( \xi _ { \alpha } - \xi _ { \alpha + 1 } \right) \cdots \left( \xi _ { \alpha } - \xi _ { n _ { \ell } } \right) }

where :math:`\alpha = 0,...,n_l`, and :math:`n_l` will be 5 to 10 (in
our case).
The implemented algorithm follows this formulation strictly. In multiple
dimensions, the shapefunctions :math:`N`, are defined as the product of
the polynomials in each dimension., i.e.
:math:`N_{kl}(\xi,\eta) = \ell^n_k(\xi)\ell^n_l(\eta).`



Derivative of the Lagrange Polynomials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As seen in the formulation section, the spectral element method also
requires the derivative of the shape functions. The derivative of the
Lagrange polynomials are defined as the original Lagrange polynomial
multiplied by the Legendre polynomial. The Legendre polynomial
:math:`P^n_i` is given by following sum

.. math:: P^n_i(\xi) = \sum_{\substack{j = 0\\ j \neq i}}^{n+1} \frac{1}{\xi - \xi_j},

where :math:`\xi_i` are the collocation points, :math:`i` is the
polynomial number with respect to the collocation point and :math:`n` is
the polynomial degree. Hence, the derivative of the Lagrange polynomial
(Eq. [eq:lagrange]) is given as

.. math:: \frac{d\ell^n_i}{d\xi} = P^n_i(\xi)\ell^n_i(\xi).

The equation becomes more elaborate when writing out all the component:

.. math:: \frac{d\ell^n_i}{d\xi} =  \sum_{\substack{k = 0\\ k \neq i}}^{n+1} \left( \frac{1}{\xi - \xi_k} \prod_{\substack{j = 0\\ j \neq i}}^{n+1} \frac{\xi - \xi_j}{\xi_i - \xi_j} \right)

This cannot be implemented as strictly as the Lagrange polynomial itself
since it would include a division by :math:`0` at :math:`\xi = \xi_k`.
However, a quick look at the denominator in the first term and the
numerator on the second term shows, that we have a cancellation when
:math:`(\xi - \xi_k) =  (\xi - \xi_j)`, or shorter, when :math:`k=j`,
such that

.. math::

   \frac{d\ell^n_i}{d\xi} =  \sum_{\substack{k = 0\\ k \neq i}}^{n+1} \left[ \prod_{\substack{j = 0\\ j \neq i}}^{n+1} \left\{ \begin{array}{ll}
           \frac{1}{\xi_i - \xi_j},  & \text{if } k=j \\
           \frac{\xi - \xi_j}{\xi_i - \xi_j}, & \text{otherwise}
       \end{array} \right. \right].

This is strictly followed in the algorithm.

The Jacobian
~~~~~~~~~~~~

The Jacobian between local and global coordinates can be defined as

.. math:: 

       \frac { \partial \mathbf { x } } { \partial \xi } = \sum _ { a = 1 } ^ { n _ { a } } \frac { \partial N _ { a } } { \partial \xi } \mathbf { x } _ { a }

In practice we use the inverse Jacobian matrix for each node to 
compute the global derivative at its location using 
``src.gll_library.dN_local_2D`` and ``src.gll_library.global_derivative``.
``dN_local_2D()`` computes the local derivative matrix containing derivatives
in each dimension of each shape function and ``global_derivative``, then takes in 
the Jacobian computes its inverse and multiplies it with the local derivative 
matrix to compute the global derivative matrix.


Global Matrix Assembly
^^^^^^^^^^^^^^^^^^^^^^

Since all elements are interconnected and dependent on each other, they
need to be assembled in a linear system of equations.

-  diagonal Mass matrix, simplifies things...

-  this, that

Then, we end up with following system of equations:

.. math:: \mathbf{M} \ddot { \mathbf{u} } + \mathbf{K} \mathbf{u} = \mathbf{f}

However, :math:`\mathbf{u}` is still dependent on time. To solve this a
simple finite difference scheme can be applied or other time marching
methods.




Time Stepping
-------------

After reduction, we get the following ordinary differential equations we
are going to solve

.. math:: M\ddot{ U } +KU=F(t)

where M is the global mass matrix and it is time invarient, the K is the
global stiff matrix and it is also time invarient. In our case, the mass
matrix M is diagonal.

Generally, there are two approaches. One way is to transform the second
order equations into first order equations by :math:`P=\dot { U }`,
then we need to solve the following equations

.. math:: \begin{pmatrix} I & O \\ O & M \end{pmatrix}\left( \begin{array}{c} \dot { U }  \\ \dot { P }  \end{array} \right) +\begin{pmatrix} O & -I \\ K & O \end{pmatrix}\left( \begin{array}{c} U \\ P \end{array} \right) =\left( \begin{array}{c} O \\ F \end{array} \right)

using. The new global matrix is still diagonal. So we can use some
explicit time schemes like explicit Euler, fourth order Runge-Kutta to
solve the equations.

The other way is to solve the second order equations directly. There are
schemes like Newmark-\ :math:`\beta` and HHT-\ :math:`\alpha`. In order
to fully exploit the advantage of diagonal matrix, we only use a special
case of Newmark-\ :math:`\beta` scheme, where :math:`\beta=0`. The
algorithm is as follows:

+----+--------------------------------------------------------+
| **Algorithm 1**                                             |
+====+========================================================+
| 1: | **Procedure** :math:`\text{NEWMARK}(U_n,V_n,A_n,dt,F)` |
+----+--------------------------------------------------------+
| 2: | :math:`U_{n+1}=U_n+V_ndt+0.5A_ndt^2`                   |
+----+--------------------------------------------------------+
| 3: | :math:`A_{n+1}=M^{-1}F_{n+1}-M^{-1}KU_{n+1}`           |
+----+--------------------------------------------------------+
| 4: | :math:`V_{n+1}=(1-\gamma)A_{n}dt+\gamma A_{n+1}dt`     |
+----+--------------------------------------------------------+
| 5: | **return** :math:`U_{n+1},V_{n+1},A_{n+1}`             |
+----+--------------------------------------------------------+



where :math:`V=\dot { U } ,A=\ddot { U }`. Typically, since the
dimensions of the matrix in Newmark scheme is smaller, the computation
is faster than indirect schemes. This is one reason why Newmark scheme
is popular in waveform modelling in seismology. So in our simulations,
we will mainly use Newmark scheme.


