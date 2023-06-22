CC = gcc
MCC = mpicc
CFLAGS = -Wall -g -c
LDFLAGS = -lm

all: libphysics.a

libphysics.a: physics.o oscserial.o utils.o thermoutils.o oscpara.o thermopara.o thermoserial.o
	ar ruv libphysics.a physics.o oscserial.o utils.o thermoutils.o oscpara.o thermopara.o thermoserial.o
	ranlib libphysics.a

physics.o: physics.c oscserial.h oscpara.h
	$(CC) $(CFLAGS) physics.c -o physics.o

oscserial.o: oscserial.c oscserial.h
	$(CC) $(CFLAGS) oscserial.c -o oscserial.o

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) utils.c -o utils.o

thermoutils.o: thermoutils.c thermoutils.h
	$(CC) $(CFLAGS) thermoutils.c -o thermoutils.o

oscpara.o: oscpara.c oscpara.h thermopara.h thermoserial.h
	$(MCC) $(CFLAGS) oscpara.c -o oscpara.o

thermopara.o: thermopara.c thermopara.h
	$(MCC) $(CFLAGS) thermopara.c -o thermopara.o

thermoserial.o: thermoserial.c thermoserial.h
	$(MCC) $(CFLAGS) thermoserial.c -o thermoserial.o

clean:
	rm -f *.o *.a

.PHONY: all clean
