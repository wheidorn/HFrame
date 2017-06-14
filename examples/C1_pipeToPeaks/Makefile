CPP=g++
LD=g++
SPECIALFLAGS=-O3 -std=c++11 -stdlib=libc++
#SPECIALFLAGS=-g -O0
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

#HFramesDir=/users/jieyu/Workarea/HFramesNewTest
HFramesDir=/Users/jieyu/Workarea/StaveQC/PkgFrameTH2/HFrames

CFLAGS = $(SPECIALFLAGS) -I. -I$(HFramesDir)/include
LFLAGS = -L. -L$(HFramesDir)/bin -lhframe

DYLD_LIBRARY_PATH=$(HFramesDir)/bin:$(DYLD_LIBRARY_PATH)

RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS) -lSpectrum 

SRC=d.cpp tool.cpp

all: $(SRC:.cpp=.o)
	$(LD) $(SRC:.cpp=.o) $(RLXX) -lMinuit -o dana

d.o: d.cpp tool.h
	$(CPP) $(RCXX) -c d.cpp

tool.o: tool.cpp tool.h
	$(CPP) $(RCXX) -c tool.cpp

clean:
	rm -f *~ *.o
	rm -f dana
