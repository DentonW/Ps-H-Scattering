FC = ifort
CC = gcc
#FFLAGS = -O2 -132 -parallel -openmp -gen-interfaces -warn interfaces
FFLAGS = -O2 -132 -parallel -openmp
LDLIBS = -llapack -lblas -lgsl -lgslcblas -lstdc++
#OBJS = PsHMain.o Minimize.o gamma.o HylleraasIntegral.o HylleraasIntegralRecursion.o dilog.o ydplog.o ygam.o ygam1s.o idx.o invsg.o leq1s.o
OBJS = PsHMain.o gamma.o HylleraasIntegral.o HylleraasIntegralRecursion.o dilog.o ydplog.o ygam.o ygam1s.o idx.o invsg.o leq1s.o DrakeYan97.o 3j.o 6j.o Asymptotic.o

PsHMain: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS)

PsHMain.o: PsHMain.f90
	$(FC) -c $(FFLAGS) PsHMain.f90

HylleraasIntegral.o: HylleraasIntegral.f90
	$(FC) -c $(FFLAGS) HylleraasIntegral.f90

HylleraasIntegralRecursion.o: HylleraasIntegralRecursion.f90
	$(FC) -c $(FFLAGS) HylleraasIntegralRecursion.f90
	
Minimize.o: Minimize.c
	$(CC) -c $(CFLAGS) Minimize.c

gamma.o: Pachucki/gamma.f
	$(FC) -c $(FFLAGS) Pachucki/gamma.f

dilog.o: Pachucki/dilog.f90
	$(FC) -c $(FFLAGS) Pachucki/dilog.f90

ydplog.o: Pachucki/ydplog.f90
	$(FC) -c $(FFLAGS) Pachucki/ydplog.f90

ygam.o: Pachucki/ygam.f90
	$(FC) -c $(FFLAGS) Pachucki/ygam.f90

ygam1s.o: Pachucki/ygam1s.f90
	$(FC) -c $(FFLAGS) Pachucki/ygam1s.f90

idx.o: Pachucki/idx.f
	$(FC) -c $(FFLAGS) Pachucki/idx.f

invsg.o: Pachucki/invsg.f90
	$(FC) -c $(FFLAGS) Pachucki/invsg.f90

leq1s.o: Pachucki/leq1s.f90
	$(FC) -c $(FFLAGS) Pachucki/leq1s.f90
	
3j.o: GeneralAsymptotic/3j.f90
	$(FC) -c $(FFLAGS) GeneralAsymptotic/3j.f90

6j.o: GeneralAsymptotic/6j.f90
	$(FC) -c $(FFLAGS) GeneralAsymptotic/6j.f90

Asymptotic.o: GeneralAsymptotic/Asymptotic.f90
	$(FC) -c $(FFLAGS) GeneralAsymptotic/Asymptotic.f90

DrakeYan97.o: GeneralAsymptotic/DrakeYan97.f90
	$(FC) -c $(FFLAGS) GeneralAsymptotic/DrakeYan97.f90

clean:
	rm -f myprog *.o *.log *.fln *~
