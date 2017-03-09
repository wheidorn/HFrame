CC=g++ # This is the main compiler
SPECIALFLAGS=-O3 -std=c++11 
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

SRCDIR=src
BUILDDIR=build
TARGETDIR=bin
TARGET=$(TARGETDIR)/libhframe.so
 
SRCEXT=cpp
SOURCES=$(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS=-g $(SPECIALFLAGS) $(ROOTCFLAGS) # -Wall 
LIB=$(ROOTLIBS) # "-L lib" # uncomment it if lib/ is not empty!
INC=-I include

$(TARGET): $(OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CC) -shared $^ -o $(TARGET) $(LIB); $(CC) -shared $^ -c -o $(TARGET) $(LIB)
#	$(CC) $^ -c -o $(TARGET) $(LIB); $(CC) $^ -c -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(INC) -fPIC -c -o $@ $<; $(CC) $(CFLAGS) $(INC) -fPIC -c -o $@ $<

clean:
	$(RM) -r $(BUILDDIR) $(TARGET); $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
