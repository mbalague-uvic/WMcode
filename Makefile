CC      = g++
CFLAGS  = -O -ffast-math -W -Wall -pedantic -ansi -Winline
LFLAGS =

mf: brain.o pool.o mfgus10Sb.o; $(CC) $(LFLAGS) brain.o pool.o mfgus10Sb.o -o mf

mfgus10Sb.o: brain.h pool.h mfgus10Sb.cpp; $(CC) $(CFLAGS) -c mfgus10Sb.cpp

brain.o: brain.h brain.cpp pool.h; $(CC) $(CFLAGS) -c brain.cpp 

pool.o: pool.h pool.cpp; $(CC) $(CFLAGS) -c pool.cpp


.PHONY: clean
clean:
	rm -f *.o
	rm -f mf
