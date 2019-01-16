# TODO


## Parallelization

## 2D --> 3D

Places in the code that need to be modified to create a generic code:

1. mesh_spec.py
    1. mesh_interp2d
        1. 277: Node numbering has to be hardcoded, since it’s defined 
                by the external mesh or better it should be an input.
        2. 288: Extension in Dimension (chi)
        3. 305: extent array in first dimension
        4. 312: extent lagrange polynomial to 3D
        5. 328: glob_gll must be of X, Y and Z
        6. 345: xbool needs to have a 3rd entry
2. gll_library.py
    1. 030: lagrange(i, x, xi)
    2. NEW: Necessary: Lagrange 3D that takes in all dimensions
    3. NEW: Necessary: Lagrange1st 3D that takes in all dimensions
    4. 153: LagrangeDerMat3D  needed
    5. 240: dN_local_3D needed
    6. 441: flattened weights change to 3D

3. Stiffness 3D dimension

4. Mass Matrix: extension to 3D?


### Independent of dimension
* Time-Schemes only take in 2D matrices and the matrices are agnostic of dimensions
* Computation of stiffness matrix and mass matrix is (more or less) independent of dimension



