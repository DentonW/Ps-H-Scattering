FC = ifort
CC = gcc
OBJDIR = obj
FFLAGS = -O2 -132 -parallel -openmp -I $(OBJDIR) -module $(OBJDIR)
LDLIBS = -llapack -lblas -lgsl -lgslcblas -lstdc++

OBJS = $(addprefix $(OBJDIR)/, PsHMain.o lu.o gamma.o HylleraasIntegral.o HylleraasIntegralRecursion.o dilog.o ydplog.o ygam.o ygam1s.o idx.o invsg.o leq1s.o DrakeYan97.o 3j.o 6j.o Asymptotic.o)

PsHMain: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: Pachucki/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: Pachucki/%.f
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: GeneralAsymptotic/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<


# http://www.gnu.org/software/make/manual/make.html#Prerequisite-Types
$(OBJS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)


clean:
	rm -f PsHMain $(OBJDIR)/*.o
