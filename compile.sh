#!/bin/sh
ftn lapack.f90 isotopomer_coordinates_conversion.f90 conversion_test.f90 -llapack
