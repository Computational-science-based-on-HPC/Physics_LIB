CC = mpicc
AR = ar
CFLAGS = -Wall -Wextra -Iinclude -lm
LDFLAGS = -L. -lmylib -lmpi

SRCS = src/physics.c src/oscserial.c src/oscpara.c src/utils.c src/thermoserial.c src/thermopara.c src/thermoutils.c
OBJS = $(SRCS:.c=.o)

LIB_NAME = myphysicslib.a

.PHONY: all clean

all: $(LIB_NAME)
%.o: %.c
    $(CC) $(CFLAGS) -c $< -o $@
    
$(LIB_NAME): $(OBJS)
    $(AR) rcs $(LIB_NAME) $(OBJS)


clean:
    rm -f $(OBJS) $(LIB_NAME)

