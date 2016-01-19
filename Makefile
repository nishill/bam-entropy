# Makefile

MKFILE     = Makefile
CC         = g++ -g -O0 -Wall -Wextra -std=gnu++11
IDIR       = ./bamtools/include/
LDIR       = ./bamtools/lib/
CFLAGS     = -I $(IDIR) 
LFLAGS     = -L $(LDIR) -lbamtools
CPPSOURCE  = main.cpp BamAlignmentIterator.cpp CompressionIterator.cpp
CPPHEADER  = BamAlignmentIterator.h CompressionIterator.h
OTHERS     = ${MKFILE} README
OBJECTS    = $(CPPSOURCE:.cpp=.o)
ALLSOURCES = ${CPPHEADER} ${CPPSOURCE} ${OTHERS}
EXECUTABLE = compress 

all : ${ALLSOURCES} ${EXECUTABLE}

%.o : %.cpp
	${CC} $(CFLAGS) $(LFLAGS) -c $<

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(CFLAGS) $(LFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(LFLAGS) $< -o $@

clean:
	- rm $(OBJECTS) ${EXECUTABLE}
