CC = gcc
MCC = mpicc
CFLAGS = -Wall -g -c
LDFLAGS = -lm
OMPFLAG = -fopenmp
SRC_DIR = src
OBJ_DIR = make
INCLUDE_DIR = include

OBJS = physics.o oscserial.o utils.o thermoutils.o oscpara.o thermopara.o thermoserial.o
all: libphysics.a

libphysics.a: $(OBJS)
	ar ruv libphysics.a $(OBJS)
	ranlib libphysics.a

physics.o: ../$(SRC_DIR)/physics.c  ../$(INCLUDE_DIR)/physics.h
	$(MCC) $(CFLAGS)  ../$(SRC_DIR)/physics.c -o  ../$(OBJ_DIR)/physics.o

oscserial.o:  ../$(SRC_DIR)/oscserial.c  ../$(INCLUDE_DIR)/oscserial.h
	$(CC) $(CFLAGS)  ../$(SRC_DIR)/oscserial.c -o  ../$(OBJ_DIR)/oscserial.o

utils.o:  ../$(SRC_DIR)/utils.c  ../$(INCLUDE_DIR)/utils.h
	$(CC) $(CFLAGS)  ../$(SRC_DIR)/utils.c -o  ../$(OBJ_DIR)/utils.o

thermoutils.o:  ../$(SRC_DIR)/thermoutils.c  ../$(INCLUDE_DIR)/thermoutils.h
	$(CC) $(CFLAGS)  ../$(SRC_DIR)/thermoutils.c -o  ../$(OBJ_DIR)/thermoutils.o

oscpara.o:  ../$(SRC_DIR)/oscpara.c  ../$(INCLUDE_DIR)/oscpara.h
	$(MCC) $(CFLAGS)  ../$(SRC_DIR)/oscpara.c -o  ../$(OBJ_DIR)/oscpara.o $(OMPFLAG)

thermopara.o:  ../$(SRC_DIR)/thermopara.c  ../$(INCLUDE_DIR)/thermopara.h
	$(MCC) $(CFLAGS)  ../$(SRC_DIR)/thermopara.c -o  ../$(OBJ_DIR)/thermopara.o $(OMPFLAG)
thermoserial.o:  ../$(SRC_DIR)/thermoserial.c  ../$(INCLUDE_DIR)/thermoserial.h
	$(CC) $(CFLAGS)  ../$(SRC_DIR)/thermoserial.c -o  ../$(OBJ_DIR)/thermoserial.o

clean:
	rm -f *.o

.PHONY: all clean
