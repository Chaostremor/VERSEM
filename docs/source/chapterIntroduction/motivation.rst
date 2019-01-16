Motivation
==========

Over the last 20 years, the seismology division of the Geoscience
Department at Princeton University has been developing a software suite
to solve the elastic wave equation using the Spectral Element method
[Komatitsch,1999]. Since the spectral element method
cannot only be used to solve the wave equation, but partial differential
equations in general, the group has started to develop solvers for other
equations as well, e.g. Maxwell’s equation for magnetic anomalies. Since
the original code was setup to only compute the elastic wave equation,
it is structured to be 100% efficient and not modular. Its organisation
is too hard to grasp, making it hard to learn, understand and especially
modify. We want to overcome this issue by using the ideas of the old
code base and create a new version of the code that is structured well
and easy to modify in terms of adding a PDE solver. We will tackle this
problem by first using the example of the wave equation. However, it will 
be taken into account that the code is supposed to not only solve for the 
wave equation but also other PDE’s.


Prior Work
==========

The work that has been done until now (referring to the existing code)
is all written in Fortran90. Since our goal is to make the code human
readable as well as modular, we will write the entire code from scratch. 
That does not mean that we want to replace the existing code, but that 
we will use it a reference for a new version. 



Bibliography
============

 * Komatitsch, Dimitri and Tromp, Jeroen. Introduction to the Spectral Element Method for Three-Dimensional Seismic Wave Propagation (1999). Geophysical Journal International,Vol. 139, Issue 3, p. 806-822, DOI: 10.1046/j.1365-246x.1999.00967. https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-246x.1999.00967.x.