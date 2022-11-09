#!/bin/sh
ftn lapack.f90 isotopologue_coordinates_conversion.f90 conversion_test.f90 -llapack
