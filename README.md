# FORTRESS

FORTRESS software package consists of a set of three main FORTRAN 90/95 codes with OpenMP parallelization and 
can be used to solve the coupled Gross-Pitaevskii equations (CGPEs) for a spin-orbit-coupled spin-1 or a spin-2 
Bose-Einstein condensate (BEC). The codes are apt to compute the static properties like the ground state energy, density
profiles, chemical potentials, and rms sizes, etc. as well as the real-time dynamics.

 
## What the user need to have on the system to use the package?

The codes use FFTW software library which is freely available on [fftw.org](http://www.fftw.org/). The software library needs to be installed 
on the system if not already installed. Please refer to [FFTW Installation and Customization](http://www.fftw.org/fftw2_doc/fftw_6.html) for the installation 
instructions. Please take note that while configuring FFTW suitable flags (--with-openmp --enable-threads) need to be
invoked to enable FFTW threading.

The system should have a FORTRAN compiler installed on it. It can be any of freely available FORTRAN compilers like
GNU's gfortran, Fortran 95 compiler f95 or to have a better performance it can be Intel's FORTRAN compiler. Testing of
the codes has been done with all these three compilers.

## Directory structure of FORTRESS

The FORTRESS package contains a Makefile and following three sub-folders:
- SRC 
- INPUT  
- OUTPUT 

### SRC has the following codes:

- cgpe1d.f90
- cgpe2d.f90
- cgpe3d.f90
- diff.f90
- simpson.f90

**cgpe3d.f90** is the main program file to solve CGPES in imaginary or realtime for a three-dimensional 
(3D) SO-coupled spin-f BEC, where f = 1 or 2. Similarly **cgpe2d.f90** and **cgpe1d.f90** are main 
program files to solve the CGPEs for quasi-two-dimensional (q2D) and quasi-one-dimensional (q1D) 
BECs, respectively. The **diff.f90** computes the first order derivative of one-dimensional function 
using nine-point Richardson's extrapolation formula, and **simpson.f90** integrates the function 
by Simpson's (1/3) rule.   

### INPUT folder has three input files: 

- input1D.dat 
- input2D.dat 
- input3D.dat

Input parameters must be provided by user for q1D, q2D, and 3D systems via these files. For realtime simulation, the
package requires an initial solution written in file INPUT/initial\_sol\_spinSPIN\_ID.dat, where SPIN\_ID = 1 for spin-1
and 2 for spin-2 systems, respectively. The contents of this file need to be in the same pattern as that of 
'solution\_file\_im\_spinSPIN\_ID.dat' as described below.

#### Description of input parameters provided in input3D.dat file

>    SPIN                   (1 for spin-1 or 2 for spin-2)  

>    NITER                  (Total number of time iterations)

>    STP                    (Number of iterations after which energy, chemical potentials, and rms sizes corresponding to each component are calculated)

>    NSTP                   (Number of iterations after which component densities and their corresponding phases are written)

>    NUMBER OF THREADS      (Number of OpenMP and FFTW3 threads)

>    NX, NY, NZ             (Number of spatial-grid along x, y and z directions, respectively; for q1D only NX and for q2D system only NX and NY are to be provided)

>    DX, DY, DZ, DT         (Spatial and temporal step-sizes; for q1D only DX and DT and for q2D DX, DY, and DT needs to provided)

>    MASS OF ATOM           (Mass of atom in AMU)

>    A0, A2, A4             (Scattering lengths corresponding to total spin
>                            channels 0, 2 and 4, respectively; A4 needs to be provided for SPIN = 2 only)

>    NATOMS                 (Total number of atoms)

>    NUX, NUY, NUZ          (Trapping frequencies in Hz along x, y, and z axes, respectively)

>    TOL                    (Tolerance)

>    SWITCH\_IM             (It is set to 1 or 0 to choose imaginary or real-time propagation)

>    GAMMAX, GAMMAY, GAMMAZ (SO-coupling strenghts along x, y, and z directions)

>    OPTION\_FPC            (1 for ferro, 2 for polar, and 3 for cyclic and will be read in the absence of SO coupling)

>    MAG                    (Magnetization which will be read only in imaginary time-propagation in the absence of SO coupling)

### OUTPUT folder contains the output data files. 

More description of the contents of this folder follows later in this file.

## Use of Makefile to create the executable

-  Compiling the cgpe1d.f90 program with the Intel Fortran compiler:
   ```sh
   make cgpe.x COMPILER=ifort DIMENSION=1D
   ```
   
   > Other options for COMPILER are gfortran, gfortran\_mkl, and f95 
   
   > Other options for DIMENSION are 2D and 3D
   
   > make cgpe.x is equivalent to COMPILER=ifort DIMENSION=3D
   
   > The COMPILER and DIMENSION variables can also be set by editing the Makefile

- Running the compiled program:
  ```sh
  ./cgpe.x
  ```
- To clean the executable, err\_msg.dat, .mod, and all object data files
  ```sh
  make clean
  ```
- To clean the executable, err\_msg.dat, .mod, all the object and output data files
  ```sh
  make cleanall
  ```

### Contents of the OUTPUT folder

After running the executable in FORTRESS folder the output files written in OUTPUT folder for SWITCH\_IM = 1 are:
 - file1\_im\_spinSPIN\_ID.dat 
 - file2\_im\_spinSPIN\_ID.dat 
 - file3\_im\_spinSPIN\_ID.dat
 - tmp\_solution\_file\_spinSPIN\_ID.dat for DIMENSION=1D and 2D
 - tmp\_solution\_file\_spinSPIN\_ID.dat, tmp\_solution\_file\_spinSPIN\_ID\_xy.dat, 
tmp\_solution\_file\_spinSPIN\_ID\_xz.dat for DIMENSION=3D 
 - solution\_file\_im\_spinSPIN\_ID.dat 
 - convergence.dat

If SWITCH\_IM = 0, 'im' in the name of some of the above file gets replaced with 're'.
 

#### Description of contents of 'file1\_im\_spinSPIN\_ID.dat' or  'file1\_re\_spinSPIN\_ID.dat'

> Various input parameters are written at the top of file

> Total norm, energy, chemical potentials of the component wave-functions and absolute values of component
wave-functions at the origin corresponding to initial or guess solution, for the transient solution obtained after 
NSTP time iterations and for the converged solution (in realtime code converged solution will simply correspond
to the solution after NITER iterations).
    

#### Description of contents of 'file2\_im\_spinSPIN\_ID.dat' or 'file2\_re\_spinSPIN\_ID.dat'

Following data is written after every STP iterations in this file.

> Time is written in first column, energy in second column, and rms sizes of the component wavefunctions in 
> the last three or five columns.



#### Description of contents of 'file3\_im\_spinSPIN\_ID.dat' or 'file3\_re\_spinSPIN\_ID.dat'.

Following data is written after every STP iterations in this file.

> Time is written in first column, norm of the individual components in subsequent three or five columns, 
> followed by the total norm and magnetization in the last two columns.

#### Description of contents of 'tmp\_solution\_file\_spinSPIN\_ID.dat'.

In this file, transient solutions is written after every NSTP iterations.

> For DIMENSION = 1D, first column has x variable, density of first component followed by its phase in next two 
> columns, respectively, then the same for rest of the components.

> For DIMENSION = 2D, first two columns would have x and y variables, respectively, density of first component 
> followed by its phase in next two columns, respectively, then the same for rest of the components.

> For DIMENSION = 3D, density of first component followed by its phase in the first two columns, respectively, 
> then the same for the rest of the components

#### Description of the contents of 'convergence.dat'

This file is written in only imaginary time propagation, i.e., when SWITCH\_IM = 1.
Following data is written after every STP iterations in this file.
> Time is written in the first column, then wavefunction convergence parameter in second column, and energy 
> convergence parameter in the last column.



#### Description of 'solution\_file\_im\_spinSPIN\_ID.dat' or 'solution\_file\_re\_spinSPIN\_ID.dat' 

> In this file, final converged solution is written in case SWITCH\_IM = 1 or solution after NITER is written 
> in case SWITCH\_IM = 0 with the same entries as in 'tmp\_solution\_file\_spinSPIN\_ID.dat'.



