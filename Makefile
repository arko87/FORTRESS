# Makefile for FORTRESS package 
# Authors: Paramjeet Banger, Pardeep Kaur, Arko Roy, & Sandeep Gautam
#
# For compilation option (1),(3) and (4) Intel's Math Kernal Library needs to be
# installed on the system if not already installed with intel's compiler. For installation
# instructions please refer to
# https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html.

SHELL = /bin/sh
COMPILER = ifort
DIMENSION = 3D
DEBUG = 
ifeq ($(COMPILER), ifort)
   FC = ifort 
   INCLUDES =
   FFLAGS = -qopenmp -mkl=parallel
   LDFLAGS = -lfftw3_omp -lfftw3 -lm
   #DEBUG = -traceback 
   #DEBUG = -O0 -g -traceback -check all -check bounds -check uninit -fp-stack-check -ftrapuv -debug all \
   #       -gen-interfaces -warn interface
endif

ifeq ($(COMPILER), gfortran)
   FC = gfortran 
   INCLUDES = -I/usr/local/include/
   FFLAGS = -ffree-form -fopenmp
   LDFLAGS = -lfftw3_omp -lfftw3 -lm -llapack -lblas
   #DEBUG = -fbacktrace 
   #DEBUG = -g3 -fbacktrace -pedantic -std=f2008 -pedantic-errors -Wall -Waliasing -Wampersand -Wc-binding-type -Wcharacter-truncation \
   #      -Wline-truncation -Wconversion -Wextra -Wimplicit-procedure -Wintrinsics-std -Wreal-q-constant -Wsurprising -Wtabs \
   #      -Wunderflow -Wintrinsic-shadow -Wunused-dummy-argument -Wunused-parameter -Walign-commons -Wfunction-elimination \
   #      -Wrealloc-lhs -Wrealloc-lhs-all -Wcompare-reals -Wtarget-lifetime -fcheck=all -fcheck=bounds -fcheck=array-temps \
   #      -ffpe-trap=zero,overflow,invalid 
endif

ifeq ($(COMPILER), gfortran_mkl)
   FC = gfortran 
   INCLUDES = -I/usr/local/include/ -I${MKLROOT}/include
   FFLAGS = -ffree-form -fopenmp -m64 -Wl,--no-as-needed 
   LDFLAGS = -lfftw3_omp -lfftw3 -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
   #DEBUG = -fbacktrace
   #DEBUG = -g3 -fbacktrace -pedantic -std=f2008 -pedantic-errors -Wall -Waliasing -Wampersand -Wc-binding-type -Wcharacter-truncation \
   #     -Wline-truncation -Wconversion -Wextra -Wimplicit-procedure -Wintrinsics-std -Wreal-q-constant -Wsurprising -Wtabs \
   #     -Wunderflow -Wintrinsic-shadow -Wunused-dummy-argument -Wunused-parameter -Walign-commons -Wfunction-elimination \
   #     -Wrealloc-lhs -Wrealloc-lhs-all -Wcompare-reals -Wtarget-lifetime -fcheck=all -fcheck=bounds -fcheck=array-temps \
   #     -ffpe-trap=zero,overflow,invalid
endif

ifeq ($(COMPILER), f95)
   FC = f95 
   INCLUDES = -I/usr/local/include/ -I${MKLROOT}/include
   FFLAGS = -ffree-form -fopenmp -m64 -Wl,--no-as-needed
   LDFLAGS = -lfftw3_omp -lfftw3 -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
   #DEBUG = -fbacktrace
   #DEBUG = -g3 -fbacktrace -pedantic -std=f2008 -pedantic-errors -Wall -Waliasing -Wampersand -Wc-binding-type -Wcharacter-truncation \
   #     -Wline-truncation -Wconversion -Wextra -Wimplicit-procedure -Wintrinsics-std -Wreal-q-constant -Wsurprising -Wtabs \
   #     -Wunderflow -Wintrinsic-shadow -Wunused-dummy-argument -Wunused-parameter -Walign-commons -Wfunction-elimination \
   #     -Wrealloc-lhs -Wrealloc-lhs-all -Wcompare-reals -Wtarget-lifetime -fcheck=all -fcheck=bounds -fcheck=array-temps \
   #     -ffpe-trap=zero,overflow,invalid
endif

ifeq ($(DIMENSION), 1D)
   SRC = ./SRC/cgpe1d.f90 ./SRC/simpson.f90 ./SRC/diff.f90
else ifeq ($(DIMENSION), 2D)
   SRC = ./SRC/cgpe2d.f90 ./SRC/simpson.f90 ./SRC/diff.f90
else
   SRC = ./SRC/cgpe3d.f90 ./SRC/simpson.f90 ./SRC/diff.f90
endif

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(INCLUDES) $(DEBUG) $(FFLAGS) -o $@ -c $<

cgpe.x: $(OBJ)
	$(FC) $(INCLUDES) $(DEBUG) $(FFLAGS) -o $@ $(OBJ) $(LDFLAGS) 

clean: 
	@rm -f .*.x err_msg.dat *.mod ./SRC/*.o 

cleanall:
	@rm -f *.x err_msg.dat *.mod ./SRC/*.o ./OUTPUT/file*dat \
         ./OUTPUT/tmp*dat ./OUTPUT/sol*dat ./OUTPUT/conv*dat
