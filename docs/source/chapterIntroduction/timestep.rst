Timestep
-----------

Timestep is done by ``class`` Tscheme which is written in `src/tscheme.py`. It can solve ODEs
with the following type:

`First order ODEs`


First order ODEs with constant mass matrix :math:`M` and constant stiffness matrix :math:`K`
and time varient force term :math:`F(t)` can be solved. Equation then has the form
:math:`M\dot{X}+KX=F(t)`.


`Second order ODEs`


First order ODEs with constant mass matrix :math:`M` ,
constant absorbing matrix :math:`C`and constant stiffness matrix :math:`K`
and time varient force term :math:`F(t)` can be solved. Equation then has the form
:math:`M\ddot{X}+C\dot{X}+KX=F(t)`. If you are using schemes for first order equations like
`euler_explicit`, `rk4`, the equations will be turned into the following first order ODEs.




.. math:: \begin{pmatrix} I & O \\ O & M \end{pmatrix}\left( \begin{array}{c} \dot { X }  \\ \dot { P }  \end{array} \right) +\begin{pmatrix} O & -I \\ K & C \end{pmatrix}\left( \begin{array}{c} X \\ P \end{array} \right) =\left( \begin{array}{c} O \\ F \end{array} \right)


If you are using newmark, then the absorbing matrix :math:`C` has to be zero, please set `C=None` when
you use the solver.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Usage`

First Specify the following arguments

- ``M:`` mass matrix, ``numpy`` [N]x[N] array

- ``C:`` absorbing matrix, default None

- ``K:`` stiffness matrix, ``numpy`` [N]x[N] array

- ``f:`` force term f(t), given t return the force vector of type ``numpy`` [N] array.

- ``t:`` time vector for stepping, t[0] should be the start time

- ``outdir:`` result output directory

- ``gamma:`` parameter for Newmark scheme, default 0.5

- ``solver:`` name of the main solver, default euler_explicit

- ``solver2:`` solver for the first few steps of the multistep method,
                    default rk4

- ``interval:`` interval for output the result, should be an integer,
                     default 1


Then create a `Tscheme` object by


``tstep=tscheme.Tscheme(M=M,K=K,f=f,t=t,outdir=outdir,solver=solver)``

Then run the solver by

``tstep.process()``



^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Details of implementation`

For all the kinds of solvers, the Tscheme will automatically turn the equations into
standard form. That is using transformation :math:`U=MX` to make mass matrix :math:`M`
unit. The `src/tsolver.py` contains all the solver thant can solve equations with
unit mass matrix.

For first order solvers, like `euler_explicit`, `rk4`, the solvers in `src/tsolver.py` can
solve equations with form :math:`\dot{U}=F(U,t)`, the equation we are going to solve is
:math:`\dot{U}+K{ M }^{ -1 }U=F(t)`, here our :math:`F(U,t)=F(t)-K{ M }^{ -1 }U`.

For second order solvers newmark, the algorithm is as follows:

+----+--------------------------------------------------------+
| **Algorithm 1**                                             |
+====+========================================================+
| 1: | **Procedure** :math:`\text{NEWMARK}(U_n,V_n,A_n,dt,F)` |
+----+--------------------------------------------------------+
| 2: | :math:`U_{n+1}=U_n+V_ndt+0.5A_ndt^2`                   |
+----+--------------------------------------------------------+
| 3: | :math:`A_{n+1}=F_{n+1}-KM^{-1}U_{n+1}`                 |
+----+--------------------------------------------------------+
| 4: | :math:`V_{n+1}=(1-\gamma)A_{n}dt+\gamma A_{n+1}dt`     |
+----+--------------------------------------------------------+
| 5: | **return** :math:`U_{n+1},V_{n+1},A_{n+1}`             |
+----+--------------------------------------------------------+

After solving :math:`U`, use transformation :math:`X={M}^{-1}U` to get :math:`X`


