#!/usr/bin/python3
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--package-prefix-hash=/home/chenzili/petsc-hash-pkgs',
    '--with-blaslapack-dir=/opt/intel/oneapi/mkl/2021.3.0',
    '--with-cc=mpiicc',
    '--with-cxx=mpiicpc',
    '--with-fc=mpiifort',
    '--with-mkl_cpardiso-dir=/opt/intel/oneapi/mkl/2021.3.0',
    '--with-mkl_pardiso-dir=/opt/intel/oneapi/mkl/2021.3.0',
    '--with-mpiexec=mpiexec.hydra',
    '--with-precision=double',
    'COPTFLAGS=-g -O',
    'CXXOPTFLAGS=-g -O',
    'FOPTFLAGS=-g -O',
    'PETSC_ARCH=arch-ci-linux-intel-mkl-single',
  ]
  configure.petsc_configure(configure_options)
