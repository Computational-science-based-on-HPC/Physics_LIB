CC = gcc
MCC = mpicc
CFLAGS = -Wall -g -c
LDFLAGS = -lm
OMPFLAG = -fopenmp
OPJS = physics.o oscserial.o utils.o thermoutils.o oscpara.o thermopara.o thermoserial.o
all: libphysics.a

libphysics.a: $(OBJS)
	ar ruv libphysics.a $(OBJS)
	ranlib libphysics.a

physics.o: physics.c physics.h
	$(CC) $(CFLAGS) physics.c -o physics.o

oscserial.o: oscserial.c oscserial.h
	$(CC) $(CFLAGS) oscserial.c -o oscserial.o

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) utils.c -o utils.o

thermoutils.o: thermoutils.c thermoutils.h
	$(CC) $(CFLAGS) thermoutils.c -o thermoutils.o

oscpara.o: oscpara.c oscpara.h
	$(MCC) $(CFLAGS) oscpara.c -o oscpara.o $(OMPFLAG)

thermopara.o: thermopara.c thermopara.h
	$(MCC) $(CFLAGS) thermopara.c -o thermopara.o $(OMPFLAG)

thermoserial.o: thermoserial.c thermoserial.h
	$(CC) $(CFLAGS) thermoserial.c -o thermoserial.o

clean:
	rm -f *.o

.PHONY: all clean
