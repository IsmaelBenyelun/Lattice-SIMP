# Lattice-SIMP

Topology Optimization algorithm for (linear elastic) truss lattice structures based on the SIMP approach

Execute using julia `main.jl` or `main-3d.jl`, or run in the command line:

```cli
$ julia main.jl
```

The script `main.jl` works as follows:

- Creates a truss lattice (unit square with two diagonal bars) given the number of joints in height `nodes_h` and the number of joints in width `nodes_w`

- Choose one of the two load cases defined, in the 2D case: `"uniaxial"` and `"cantilever"`

- Performs a compliance minimization problem subjected to maximum volume `V_max` (or `volfrac`) using MMA (from `NLOpt`)

- Penalization parameter is selected via the variable `p`, and a sensitivity filter is of a radius `R' cone kernel is imposed

- Paraview output files (`.vtu`) are generated in an `output-files` folder (which is automatically made)

The script `main-3d.jl` is a 3d generalization of the previous case. Coordinates and connectivities are required to be in the corresponding `.csv` files in the `data` folder. The load case for the current optimization is hard-coded (TBD)
