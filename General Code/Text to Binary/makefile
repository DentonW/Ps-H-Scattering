FC = g++
FFLAGS = -O2
LDLIBS = 
OBJS = Binary.o

Binary: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDLIBS)

Binary.o: Binary.cpp
	$(FC) -c $(FFLAGS) Binary.cpp

	
clean:
	rm -f myprog *.o *.log *.fln *~
