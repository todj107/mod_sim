CC = gcc

CFLAGS = -g -O3		# debugging information, optimization level
CPPFLAGS = -I.		# Flags to the C preprocessor, here include files from the current directory
LOADLIBES = -lm		# Libraries to load, here libm - the math library
VPATH = ../src		# list of directories that make should search


# List of object files
OBJS = sim.o ran.o common.o config.o measure.o dynamics.o

# Implicit rule - link together object files to an executable file
sim: ${OBJS}

# The object files also depend on some header files and Makefile
${OBJS}: Makefile sim.h define.h ran.h
