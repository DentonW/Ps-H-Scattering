FC = icpc
FFLAGS = -O3 -I/usr/include/i386-linux-gnu -I/usr/include/openmpi -xHost -openmp -I/opt/intel/mkl/include
#LDLIBS = -lmkl_intel_thread -lmkl_core -lmkl_lapack95 -lmkl_sequential -lm -lmkl_intel -lmkl_intel_thread -lmkl_blas95
LDLIBS = -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -lmkl_sequential -lm -lmkl_intel_lp64 -lboost_system -lboost_filesystem
OBJS = Phase\ Shift.o

Phase: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS) -L$MKLROOT/lib/em64t -L/opt/intel/mkl/lib/intel64
	
	#$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS) -L$MKLROOT/lib/ia32 -L/opt/intel/mkl/lib/ia32

Phase\ Shift.o: Phase\ Shift.cpp
	$(FC) -c $(FFLAGS) Phase\ Shift.cpp
	
clean:
	rm -f Phase *.o
