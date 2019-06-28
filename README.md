# VMECTools

## vmec2pest
Transforms a VMEC equilibrium into PEST coordinates.  The magnetic equilibrium is specified in the PEST coordinate system $\mathbf{B} = \nabla x^1 \times \nabla x^2$ where $x^1 = \psi_\text{toroidal}/\psi_\text{LCFS}$ is the normalized toroidal flux and $x^2 = \alpha = \that -\iota \zeta$ is the magnetic field line label parameterized by the straight-field-line coordinates $\theta$ and $\zeta$.  For the PEST system, $\zeta$ is chosen to be the geometric toroidal angle. 


### Building
Building `vmec2pest` requires the following packages
  - Fortran compiler with Fortran 2008 support
  - [NetCDF libraries](https://www.unidata.ucar.edu/software/netcdf/)
  - An optimized BLAS/LAPACK installation (such as [OpenBLAS](https://github.com/xianyi/OpenBLAS))

The `NETCDF_F_DIR`, `NETCDF_C_DIR` and `BLAS_DIR` directories need to be set to the appropriate path.  To make `vmec2pest`, in the top-level directory, execute
```
make v2p
```
which will then build the `mini_libstell` library if not built and the `vmec2pest` executable. 

### Usage
Copy the template input file from the `input/` directory to the top-level directory.  The following lines in the input file need to be set.

```
geomdir = 'absolute/path/to/VMEC/wout.nc/file'
surfaces = 0.1, 0.2, ... Array of surfaces on which to compute the transformation
n_field_lines = # of field lines on which to compute the transformation
n_parallel_pts = # of points in the parallel direction on which to compute the transformation
```

`vmec2pest` produces by default a `pest` file for each specified surface.  The file contains the metric coefficients, magnetic field strength, Jacobian, curavature drift components and parallel derivative of B as a function of field line label and field-line-following coordinate.

To execute, in the top-level directory, simply run
```
./vmec2pest
```
If the `outdir` directory is not specified in `vmec2pest.inp`, the files will be written to the top-level directory. 
