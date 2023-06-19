cmake_minimum_required(VERSION 3.25)
project(physics C)

set(CMAKE_C_STANDARD 11)

add_library(physics physics.c oscserial.c oscserial.h oscpara.h oscpara.h oscpara.c utils.c utils.h examples/main.c thermoserial.c thermoserial.h thermoserial.h thermopara.c thermopara.h)
