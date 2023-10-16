CC      = g++
CFLAGS  = -O -ffast-math -W -Wall -pedantic -ansi -Winline
LFLAGS =

mf: brain.o pool.o mfmain.o; $(CC) $(LFLAGS) brain.o pool.o mfmain.o -o mf

mfmain.o: brain.h pool.h mfmain.cpp; $(CC) $(CFLAGS) -c mfmain.cpp

brain.o: brain.h brain.cpp pool.h; $(CC) $(CFLAGS) -c brain.cpp 

pool.o: pool.h pool.cpp; $(CC) $(CFLAGS) -c pool.cpp


.PHONY: clean
clean:
	rm -f *.o
	rm -f mf
