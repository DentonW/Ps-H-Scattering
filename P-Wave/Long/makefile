FC = mpicc
FFLAGS = -O3 -I/opt/intel/mkl/include -openmp -parallel -fp-model precise -fp-model extended -cc=icpc
LDLIBS = -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lgsl -lgslcblas -lstdc++
OBJS = Ps-H\ Scattering.o Short-Range.o Long-Range.o Phase\ Shift.o Gaussian\ Integration.o Vector\ Gaussian\ Integration.o tinyxml2.o

PWaveScattering: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS) -L$MKLROOT/lib/intel64

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
	rm -f PWaveScattering *.o
