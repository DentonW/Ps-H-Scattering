FC = gfortran
FFLAGS = -fopenmp
LDLIBS = -llapack
OBJS = DrakeYan97.o Asymptotic.o 3j.o 6j.o

DrakeYan97: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS)

DrakeYan97.o: DrakeYan97.f90
	$(FC) -c $(FFLAGS) DrakeYan97.f90

Asymptotic.o: Asymptotic.f90
	$(FC) -c $(FFLAGS) Asymptotic.f90

3j.o: 3j.f90
	$(FC) -c $(FFLAGS) 3j.f90

6j.o: 6j.f90
	$(FC) -c $(FFLAGS) 6j.f90

clean:
	rm -f *.log *.fln *.o *~
