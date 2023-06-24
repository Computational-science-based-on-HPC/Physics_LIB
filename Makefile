CC = gcc
MCC = mpicc
CFLAGS = -Wall -g -c
LDFLAGS = -lm
OMPFLAG = -fopenmp
SRC_DIR = src
OBJ_DIR = obj
INCLUDE_DIR = include

OBJS = $(addprefix $(OBJ_DIR)/, physics.o oscserial.o utils.o thermoutils.o oscpara.o thermopara.o thermoserial.o)
all: libphysics.a

libphysics.a: $(OBJS)
	ar ruv libphysics.a $(OBJS)
	ranlib libphysics.a

$(OBJ_DIR)/physics.o: $(SRC_DIR)/physics.c  $(HED_DIR)/physics.h
	$(CC) $(CFLAGS)  $(SRC_DIR)/physics.c -o  $(OBJ_DIR)/physics.o

$(OBJ_DIR)/oscserial.o:  $(SRC_DIR)/oscserial.c  $(HED_DIR)/oscserial.h
	$(CC) $(CFLAGS)  $(SRC_DIR)/oscserial.c -o  $(OBJ_DIR)/oscserial.o

$(OBJ_DIR)/utils.o:  $(SRC_DIR)/utils.c  $(HED_DIR)/utils.h
	$(CC) $(CFLAGS)  $(SRC_DIR)/utils.c -o  $(OBJ_DIR)/utils.o

$(OBJ_DIR)/thermoutils.o:  $(SRC_DIR)/thermoutils.c  $(HED_DIR)/thermoutils.h
	$(CC) $(CFLAGS)  $(SRC_DIR)/thermoutils.c -o  $(OBJ_DIR)/thermoutils.o

$(OBJ_DIR)/oscpara.o:  $(SRC_DIR)/oscpara.c  $(HED_DIR)/oscpara.h
	$(MCC) $(CFLAGS)  $(SRC_DIR)/oscpara.c -o  $(OBJ_DIR)/oscpara.o $(OMPFLAG)

$(OBJ_DIR)/thermopara.o:  $(SRC_DIR)/thermopara.c  $(HED_DIR)/thermopara.h
	$(MCC) $(CFLAGS)  $(SRC_DIR)/thermopara.c -o  $(OBJ_DIR)/thermopara.o $(OMPFLAG)

$(OBJ_DIR)/thermoserial.o:  $(SRC_DIR)/thermoserial.c  $(HED_DIR)/thermoserial.h
	$(CC) $(CFLAGS)  $(SRC_DIR)/thermoserial.c -o  $(OBJ_DIR)/thermoserial.o

# clean:
# 	rm -f *.o

.PHONY: all clean
