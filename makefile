propagate: propagate.f90
# 	gcc -o hellomake hellomake.c hellofunc.c -I.
	f2py -m propagate --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c propagate.f90
