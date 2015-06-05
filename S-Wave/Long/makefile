FC = mpiCC
FFLAGS = -O2 -I/opt/intel/Compiler/11.1/069/mkl/include -I/opt/intel/Compiler/11.1/073/mkl/include -I/opt/intel/Compiler/11.1/074/mkl/include -I/usr/include/openmpi -xHost -openmp -I/opt/intel/mkl/include -I/usr/include/i386-linux-gnu -I/opt/intel/mkl/include 
#LDLIBS = -lmkl_core -lmkl_lapack95 -lmkl_sequential -lm -lmkl_intel -lmkl_intel_thread -lmkl_blas95 
LDLIBS = -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lstdc++
OBJS = Ps-H\ Scattering.o Short-Range.o Long-Range.o Phase\ Shift.o Gaussian\ Integration.o Vector\ Gaussian\ Integration.o tinyxml2.o

Scattering: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS) -L$MKLROOT/lib/em64t 
	#$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS) -L$MKLROOT/lib/em64t -L/opt/intel/mkl/lib/intel64 -L/opt/intel/Compiler/11.1/069/mkl/lib/em64t -L/opt/intel/Compiler/11.1/073/mkl/lib/em64t -L/opt/intel/Compiler/11.1/074/mkl/lib/em64t -L/opt/intel/mkl/lib/ia32

Ps-H\ Scattering.o: Ps-H\ Scattering.cpp
	$(FC) -c $(FFLAGS) Ps-H\ Scattering.cpp

Gaussian\ Integration.o: Gaussian\ Integration.cpp
	$(FC) -c $(FFLAGS) Gaussian\ Integration.cpp

Phase\ Shift.o: Phase\ Shift.cpp
	$(FC) -c $(FFLAGS) Phase\ Shift.cpp

Long-Range.o: Long-Range.cpp
	$(FC) -c $(FFLAGS) Long-Range.cpp

Short-Range.o: Short-Range.cpp
	$(FC) -c $(FFLAGS) Short-Range.cpp

Vector\ Gaussian\ Integration.o: Vector\ Gaussian\ Integration.cpp
	$(FC) -c $(FFLAGS) Vector\ Gaussian\ Integration.cpp

tinyxml2.o: tinyxml2.cpp
	$(FC) -c $(FFLAGS) tinyxml2.cpp
	
clean:
	rm -f Scattering *.o
