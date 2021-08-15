#!/usr/bin/env bash
rm a.out
ifort gaussquad.f90
./a.out
rm -r debugging
mkdir debugging

ifort Tulli13_1990.f90 progTulli12_1990.f90 -ipo -O3 -no-prec-div -static-intel -fp-model fast=2 -xHost -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

cp a.out debugging
cp fort.23 debugging
cp a.out debugging
cp parallel_script_v4.py debugging
cp raw_x.txt debugging
cp raw_w.txt debugging
cp job.sh debugging
cd debugging

python3 parallel_script_v4.py 



