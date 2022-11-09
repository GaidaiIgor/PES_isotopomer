Fortran modules for conversion between Euler coordinates of two rigid rotor molecules and Euler coordinates of their isotopomers.
Can be used to convert the Potential Energy Surface of the main isotopomer to rare isotopomers. See the comments in the source files for more detailed description.

`isotopologue_coordinates_conversion.f90` is the main conversion module. 

`constants.f90` and `lapack.f90` are auxiliary modules.

`conversion_test.f90` is an example driver for the conversion module.

Compilation is straightforward, see `compile.sh` for an example.
