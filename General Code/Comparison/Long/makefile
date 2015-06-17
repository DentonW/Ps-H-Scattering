FC = g++
FFLAGS = 
LDLIBS = -lIL -lILU
OBJS = Compare.o

Compare: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS)

Compare.o: Compare.cpp
	$(FC) -c $(FFLAGS) Compare.cpp
	

	
clean:
	rm -f Compare *.o
