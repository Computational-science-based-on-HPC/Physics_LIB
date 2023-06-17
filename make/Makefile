# need to re-check
# in SecondLib/Makefile
CC=gcc
CFLAGS= -Wall -Wextra -std=c11
BUILD_DIR = ../build
LIBRARY_NAME = myphysics
OBJ = $(BUILD_DIR)/library_source.o
LIBRARY = $(BUILD_DIR)/lib$(LIBRARY_NAME).a

## add your other stuff for make
all:  $(LIBRARY)
$(OBJ): $(SRC) ../include/physics.h ../include/oscserial.h ../include/oscpara.h
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -g -c $< -o $@
$(LIBRARY): $(OBJ)
	ar rcs $@ $<
clean:
	rm -rf $(BUILD_DIR)
.PHONY: all clean

