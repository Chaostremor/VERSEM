Overview
+++++++++++++++++++++++++++++++


Motivation
===============================

Over the last 20 years, the seismology division of the Geoscience Department at 
Princeton University has been developing a software suite to solve the elastic 
wave equation using the Spectral Element method \citep{Komatitsch1999}. 
Since the spectral element method cannot only be used to solve the wave equation, 
but partial differential equations in general, the group has started to develop 
solvers for other equations as well, e.g. Maxwell's equation for magnetic anomalies.
Since the original code was setup to only compute the elastic wave equation, it 
is structured to be 100\% efficient and not modular.
Its organisation is too hard to grasp, making it hard to learn, understand and 
especially modify.
We want to overcome this issue by using the ideas of the old code base and create
a new version of the code that is structured well and easy to modify in terms of 
adding a PDE solver. We will tackle this problem by first using the example of 
the wave equation because we have a benchmark of a working code. However, it will 
be taken into account that the code is supposed to not only solve for the wave 
equation but also other PDE's.

Prior Work
====================================

The work that has been done until now (referring to the existing code) is all written 
in Fortran90. Since our goal is to make the code human readable as modular, we will 
write the code from scratch. That does not mean that we will forget the existing code, 
but that we will use it a reference for a new version.

The Solver
====================================

The numerical solution of the wave equation using the Spectral Element Method 
(SEM), much like the Finite Element Method (FEM), is based on the weak form of 
the equation. For a discrete set of elements and corresponding nodes, PDEs can 
be solved using a linear system.



    The Weak Form of the Elastic Wave Equation
    ------------------------------------------------------------------------------
    is based on the weak solution of the elastic wave equation.
    \begin{equation}
        \rho\: \partial^2_t\mathbf{u} = \nabla \cdot \left( \mathbf{c}\, \colon\!  \nabla \mathbf{u} \right) + \mathbf{f},
    \end{equation}
    where (from left to right) $\rho$ is the density, $\partial_t$ is the differential operator with respect to time, $\mathbf{u}$ is the displacement, $\nabla$ is the gradient operator in space, $\mathbf{c}$ is the elastic tensor (Hooke's law: $\mathbf{T} = \mathbf{c}\, \colon\!  \nabla \mathbf{u} $) and $\mathbf{f}$ is an external force acting on the system. The weak form of the elastic wave equation can then be formed by the multiplying the elastic wave equation with a test function $\hat{\phi}$ and integrating it over the domain $\Omega$:
    \begin{equation}
        \int_{\Omega} \hat{\phi} \: \rho\: \partial^2_t \mathbf{u} \: dV = \int_{\Omega} \hat{\phi} \: \nabla \cdot \left( \mathbf{c}\, \colon\!  \nabla \mathbf{u} \right) \: dV + \int_{\Omega} \hat{\phi} \: \mathbf{f} \: dV
    \end{equation}
    By using Green's theorem, we can simplify the double spatial gradient:

    \begin{equation}
        \int_{\Omega} \hat{\phi} \: \rho\: \partial^2_t \mathbf{u} \: dV =  \oint_{\partial\Omega} \hat{\phi} \:\mathbf{T} \cdot \mathbf{\hat{n}} \: dS - \int_{\Omega} \nabla \hat{\phi} \cdot \left( \mathbf{c}\, \colon\!  \nabla \mathbf{u} \right) \: dV + \int_{\Omega} \hat{\phi} \: \mathbf{f} \: dV
    \end{equation}

..
    Now, if our domain can be split into elements $\mathbf{e}$ which are formed by discrete data points we can interpolate the interior of each element the using shape functions. 
    Particularly, within such elements, we can assume that the displacement $\mathbf{u}$ can be expressed in form of a sum discrete points multiplied by connecting shape functions, such that:

    \begin{equation}\label{eq:discrete_disp}
        u^i(\mathbf{x},t) =  \phi _ { 1 } ( \mathbf { x } ) u _ { 1 }^i(t) + \phi _ { 2 } ( \mathbf { x } ) u _ { 2 }^i(t) + \cdots = \sum _ { l = 1 } ^ { N n } \phi_{l}(\mathbf{x})u_{l}^i(t),
    \end{equation}
    where $Nn$ is the number of nodes per element and $i$ denotes the component of the displacement vector.
    \begin{equation}
        \sum_{l=1}^{Nn} \int_{\Omega} \hat{\phi} \: \rho\: \partial^2_t \phi_{l}(\mathbf{x})u_{l}^i(t) \: dV  + \int_{\Omega} \nabla \hat{\phi} \cdot \left( \mathbf{c}\, \colon\!  \nabla \phi_{l}(\mathbf{x})u_{l}^i(t) \right) \: dV = \oint_{\partial\Omega} \hat{\phi} \:\mathbf{T} \cdot \mathbf{\hat{n}} \: dS  + \int_{\Omega} \hat{\phi} \: \mathbf{f} \: dV  
    \end{equation}
    If we choose $\hat{\phi}(\mathbf{x})$ to be $\phi_{k}(\mathbf{x})$ for every $k$ from 1 to $Nn$, we can further reduce the problem to 

    \begin{align}
        \sum_{l=1}^{Nn} & \int_{\Omega} \phi_{k} \: \rho\: \partial^2_t \phi_{l}\: u_{l}^i(t) \: dV  \nonumber\\
        &+ \int_{\Omega} \nabla \phi_{k} \cdot \left( \mathbf{c}\, \colon\!  \nabla \phi_{l} \: u_{l}^i(t) \right) \: dV
        \nonumber\\
        &\quad -  \int_{\Omega} \phi_{k} \: \mathbf{f}  \: dV = \oint_{\partial\Omega} \phi_{k} \:\mathbf{T} \cdot \mathbf{\hat{n}} \: dS    
    \end{align}

    As given by the above equation the factors $u^i_l$ are completely independent of the integration over the domain $\Omega$. Hence, we can write

    \begin{align}
        \sum_{l=1}^{Nn} & \int_{\Omega} \rho\: \phi_{k} \:  \phi_{l} \: dV \: \partial^2_tu_{l}^i(t) \nonumber\\
        &+ \int_{\Omega} \nabla \phi_{k} \cdot \left( \mathbf{c}\, \colon\!  \nabla \phi_{l} \right) \: dV \: u_{l}^i(t)
        \nonumber\\
        &\quad - \int_{\Omega} \phi_{k} \: \mathbf{f}   \: dV = \oint_{\partial\Omega} \phi_{k} \:\mathbf{T} \cdot \mathbf{\hat{n}} \: dS  
    \end{align}

    \subsection{The Shape Functions}

    For each element, we can convert the global coordinates $\mathbf{x}$ to local coordinates that range from $( \xi , \eta ) , \quad - 1 \leq \xi \leq 1,- 1 \leq \eta \leq 1$. 
    The conversion from local to global coordinates is governed by the transformation using so called shape functions such that
    \begin{equation}
    \mathbf { x } ( \xi , \eta ) = \sum _ { a = 1 } ^ { n _ { a } } N _ { a } ( \xi , \eta ) \mathbf { x } _ { a }.
    \end{equation}
    The SEM's signature shape functions are the Lagrange Polynomials $l^{n_l}_\alpha$:
    \begin{equation}
        \ell _ { \alpha } ^ { n _ { \ell } } ( \xi ) = \frac { \left( \xi - \xi _ { 0 } \right) \cdots \left( \xi - \xi _ { \alpha - 1 } \right) \left( \xi - \xi _ { \alpha + 1 } \right) \cdots \left( \xi - \xi _ { n _ { \ell } } \right) } { \left( \xi _ { \alpha } - \xi _ { 0 } \right) \cdots \left( \xi _ { \alpha } - \xi _ { \alpha - 1 } \right) \left( \xi _ { \alpha } - \xi _ { \alpha + 1 } \right) \cdots \left( \xi _ { \alpha } - \xi _ { n _ { \ell } } \right) }
    \end{equation}
    where $\alpha = 0,...,n_l$, and $n_l$ will be 5 to 10 (in our case).

    \vspace{2ex}
    \large{These instructions will be elaborated...}

    \begin{itemize}
        \item using shape functions and the transformation with Jacobian to get from global to local integration
        \item Gauss-Lobatto-Legendre Quadrature to approximate the integrals
        \item 
    \end{itemize}

    \subsection{Global Matrix Assembly}\label{sec:gmatassembly}
    Since all elements are interconnected and dependent on each other, they need to be assembled in a linear system of equations.

    \begin{itemize}
        \item diagonal Mass matrix, simplifies things...
        \item this, that
    \end{itemize}

    Then, we end up with following system of equations:

    \begin{equation}
        \mathbf{M} \ddot { \mathbf{u} } + \mathbf{K} \mathbf{u} = \mathbf{f}
    \end{equation}

    However, $\mathbf{u}$ is still dependent on time. To solve this a simple finite difference scheme can be applied or other time marching methods.

    \subsection{Time Marching}

    The linear system of equations can then be solved using a simple finite difference scheme, where 

    \begin{equation}
    \mathbf { u } ( t + d t ) = d t ^ { 2 } \left( \mathbf { M } ^ { T } \right) ^ { - 1 } \left[ \mathbf { f } - \mathbf { K } ^ { T } \mathbf { u } \right] + 2 \mathbf { u } - \mathbf { u } ( t - d t ).
    \end{equation}

    While we will focus on this scheme first of all, we want to allow for other time schemes as well, e.g. Crank-Nicholson, to allow for more robust time marching using other PDEs.

    \section{Parallelisation}

    There are two types of parallelisation that can possibly be done. 
    The first one is the embarassingly parallel job of integration at element level. 
    This should be relatively simple implement. 
    However, the integration at element level is not the bottleneck of the computation, the bottleneck is the solving of the linear system when marching in time.
    The second type of parallelisation, which is often used in such cases, is domain decomposition and synchronisation.
    Meaning the domain is split into subdomains which overlap and are synchronised after each time-step.





    \section{Organization}

    Due to our decision to make this code as easy to understand and modular as possible, we decided to write the code in Python.
    We will however keep the option open to write certain 'under the hood functions' in C or C++, so that there is a possibility of a quick speed up.


    \subsection{Input}
    The input parameter files will be in YAML format. 
    They include information of physics and simulation of this problem. For physics, it contains the following information, solver type, mesh type, equation type, Gauss - Lobatto - Legendre Points (GLL points), post-processing software and possibly more.
    The mesh will be created externally.

    \subsubsection{Parameter File}

    A YAML parameter file will be given to provide necessary settings. Parameters will be stored as a dictionary in a configuration object, which will include methods like get(), set(), has(), etc. Paremeters include:
    \begin{enumerate}[itemsep=0em]
        \item Compute SH or P-SV wave
        \item Boundary conditions
        \item Interval of output snapshots
        \item Number of timesteps
        \item Length of one timestep
    \end{enumerate}

    \subsubsection{Source File}
    The source file is a YAML file which contains the following information:
    \begin{enumerate}[itemsep=0em]
        \item Location
        \item Frequency
        \item Source time function (either built-in SFT or points to an external STF file)
        \item Starting time
        \item Angle
        \item Moment tensor
        \item Amplification factor
    \end{enumerate}
    \subsubsection{Station File}
    The station file is a YAML file which contains the location information.

    \subsubsection{Mesh File}
    The mesh file will be created by an external meshing software.
    In this project, we will use Cubit.
    Cubit outputs an Exodus file (\texttt{.e}) which will be used to create a model object

    \section{Pre-Processing}
    The GLL points contain the total number of points and the weight at each point to approximate the integral with the weighted sum. We will use a library of hard-coded GLL points and weights for the Quadrature calculating them.  

    For each node, the Jacobian has to be calculated and saved since the locations of the nodes are not changing. 
    This is necessary to convert the coordinates from the global coordinate system to the local (elemental) coordinate system on each for the local GLL quadrature (as described in section \ref{sec:gmatassembly}).



    \subsection{\texttt{Model\_object}}
    From input parameters that described the model as well as the input mesh file are converted into a \texttt{model\_object}.
    A complete model consists of an Exodus (\texttt{.e}) file, which is produced using an external mesher and a file that describes the material. The Exodus file will be read as a \(2\times N\) array,  (\(N\) is the number of points, which contains the location of each point. The material file is a \(3\times N\) binary array (3 parameters are P-wave velocity, S-wave velocity and density). A model object that contains the arrays will be created. The \texttt{model\_object} contains methods to iterate through points and get neighbouring points ({\color{red}{I think the exodus file contains neighbouring points}}).
    

    \subsection{Output}
    \subsubsection{Seismograms}
    Seismograms records the velocity or displacement as a function of time. For each station, a binary array of velocity or displacement will be stored. The input station file will be copied to the output directory so that the output directory contains the complete seismogram information.
    \subsubsection{Wavefield Snapshots}
    Snapshots of the wavefield (either velocity or displacement) can be configured as an output. The snapshots will be stored as a \(1\times N\) binary arrays. The Exodus file of the model will be copied to the output directory, as well, so that the output contains all necessary information.

    \subsection{Post-Processing}
    The post-processing software could be Paraview, Matlab or Matplotlib in Python, which will be the default. The post-processing may not be of priority in this project.

    \subsection{Unit Testing}

    The unit test of this project will be implemented with Jenkins. Everyone will submit his own tests to the server independently, which means all the group members will do the unit test in the sections he is responsible for.

    \section{Schedule And Division of Labor}

    \subsection{Schedule}

    \subsubsection{11 Dec, 2018: Prototype}

    \subsubsection{31 Dec, 2018: Finish the main solver}

    \subsubsection{9 Jan, 2019: Final version}


    \subsection{Division of Labor}

    \textbf{Congyue Cui:} Coding of the input and output design 
    \\[8pt]
    \textbf{Chao Song:} Coding of the visualisation and details
    \\[8pt]
    \textbf{Fan Wu:} Coding of the main solver
    \\[8pt]
    \textbf{Lucas Sawade:} Coding of the main solver
    \\[8pt]
    \textbf{Srijan Bharati Das:} Coding of the main solver
    \\[8pt]
    \textbf{ALL:} Skeleton design; unit test; Documentation


    \subsection{Building the Automatic Documentation with Sphinx}



asddasda
