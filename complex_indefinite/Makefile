.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES = 
C_OBJS    = 
C_INCLUDES =
C_LIBS    = -lm
EXAMPLES = example_complex_indefinite
C_EXAMPLES_OBJ = example_complex_indefinite.o

CCX = g++
# Place your Eigen directory after the -I!
CXXFLAGS = -O2 -g -Wall -std=c++11 -I../Eigen/



#============================================================

all: $(EXAMPLES)


%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(*:=.cpp)

example_complex_indefinite :  $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)


#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(INCLUDES) $(EXAMPLES)

clean: 
	rm -f $(TARGET) $(C_OBJS) $(EXAMPLES) $(C_EXAMPLES_OBJ)

TAGS:	$(ALL_SOURCES)
	etags $(ALL_SOURCES)

tar: $(ALL_SOURCES) 
	tar cvf $(TARGET).tar $(ALL_SOURCES) 

ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


