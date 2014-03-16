# First make-file
CC=clang++
CFLAGS=-Wall #-O4 -fopenmp
IFLAGS=-I ../libpng/include
LFLAGS=-l png
LDFLAGS=-L ../libpng/lib

all:
	$(CC) *.cpp $(CFLAGS) $(LDFLAGS) $(LFLAGS) $(IFLAGS) -o ritter

debug:
	$(CC) *.cpp $(CFLAGS) $(LDFLAGS) $(LFLAGS) $(LFLAGS) -g -o ritter
