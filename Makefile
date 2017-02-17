CC=g++ # This is the main compiler
SPECIALFLAGS=-O3 -std=c++11 
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

SRCDIR=src
BUILDDIR=build
TARGET=bin/dana
 
SRCEXT=cpp
SOURCES=$(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS=-g $(SPECIALFLAGS) $(ROOTCFLAGS) # -Wall
LIB=-L lib $(ROOTLIBS)
INC=-I include

$(TARGET): $(OBJECTS)
	$(CC) $^ -c -o $(TARGET) $(LIB); $(CC) $^ -c -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	$(RM) -r $(BUILDDIR) $(TARGET); $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
