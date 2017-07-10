CC=g++
#CFLAGS=-std=c++11 -Wall -O3  -I/home/matta/armaheadlib/usr/include 
LDFLAGS=-larmadillo #-L/home/matta/armaheadlib/usr/lib 
SOURCES=main.cpp PseudoHB.cpp Sommer.cpp testing.cpp HBRandom.cpp algorithm.cpp
HEADERS=PseudoHB.h Sommer.h testing.h HBRandom.h algorithm.h latticeconfig.h
all: pseudoHB

pseudoHB: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(SOURCES) -o pseudoHB $(LDFLAGS)

clear:
	rm pseudoHB
