CC=g++ # This is the main compiler
SPECIALFLAGS=-O3 -std=c++11 -stdlib=libc++
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

SRCDIR=src
BUILDDIR=build
TARGET=bin/dana
 
SRCEXT=cpp
SOURCES=$(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS=-g $(SPECIALFLAGS) $(ROOTCFLAGS) # -Wall
#-pthread -lmongoclient  -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt
LIB=-L lib $(ROOTLIBS)
INC=-I include

#set noet in vim to type tab to avoid expandtab

$(TARGET): $(OBJECTS)
	$(CC) $^ -c -o $(TARGET) $(LIB); $(CC) $^ -c -o $(TARGET) $(LIB)
#$(CC) $^ -o $(TARGET) $(LIB); $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	$(RM) -r $(BUILDDIR) $(TARGET); $(RM) -r $(BUILDDIR) $(TARGET)

## Tests
#tester:
#  $(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester
#
## Spikes
#ticket:
#  $(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean
