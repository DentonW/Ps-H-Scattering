FC = mpif90
#FC = gfortran
#FC = ifort
#FFLAGS = -r8
#FFLAGS   = -fdefault-real-8 -O2 -fopenmp
#FFLAGS = -fopenmp -O2
#FFLAGS = -openmp -parallel -O3 -xHost
#FFLAGS = -openmp -xHost -heap-arrays -g
FFLAGS = -openmp -xHost
#LDLIBS = -lmkl_intel_thread -lmkl_lapack95 -lmkl_core -lmkl_sequential -lmkl_blacs -lmkl_blas95 -lmkl_lapack
#LDLIBS = -lmkl_core -lmkl_lapack95 -lmkl_sequential -lm -lmkl_intel -lmkl_intel_thread -lmkl_blas95
#LDLIBS = -lmkl_lapack -lmkl_intel_thread -lmkl_core -lmkl_sequential
#LDLIBS = -llapack -lblas
LDLIBS = -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential
OBJS = Energy.o invsg.o leq1s.o

Energy: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS) -L$MKLROOT/lib/em64t -L/opt/intel/Compiler/11.1/064/mkl/lib/em64t

Energy.o: Energy.f90
	$(FC) -c $(FFLAGS) Energy.f90

invsg.o: invsg.f90
	$(FC) -c $(FFLAGS) invsg.f90
	
leq1s.o: leq1s.f90
	$(FC) -c $(FFLAGS) leq1s.f90
	
	
clean:
	rm -f myprog *.o *.log *.fln *~
