# First make-file
CC=clang++
GCC=g++
CFLAGS=-Wall 
AFLAGS=-Wall -fopenmp
IFLAGS=-I ../libpng/include
AIFLAGS=-I /usr/include/libpng12
LFLAGS=-l png
LDFLAGS=-L ../libpng/lib

all:
	$(CC) *.cpp $(CFLAGS) $(LDFLAGS) $(LFLAGS) $(IFLAGS) -o ritter

debug:
	$(CC) *.cpp $(CFLAGS) $(LDFLAGS) $(LFLAGS) $(IFLAGS) -g -o ritter

aviary:
	$(GCC) *.cpp $(AFLAGS) $(LFLAGS) $(AIFLAGS) -oritter -D USE_OPEN_MP
