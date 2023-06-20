# in SecondLib/Makefile
CC=gcc
CFLAGS= -Wall -g -c
## add your other stuff for make
all: physics.o
physics.o: physics.c oscserial.h oscpara.h
	$(CC) $(CFLAGS) -c physics.c -o physics.o

oscserial.o: oscserial.c oscserial.h
	$(CC) $(CFLAGS) -c oscserial.c -o oscserial.o

oscpara.o: oscpara.c oscpara.h
	$(CC) $(CFLAGS) -c oscpara.c -o oscpara.o

libphysics.a: physics.o
	ar ruv libphysics.a physics.o
	ranlib physics.a
clean:
	rm -f *.o *.a
.PHONY: all clean