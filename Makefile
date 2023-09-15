CXX	= g++
RM	= rm -f
CSRC	= io_png.c
CXXSRC	= libdenoising.cpp lib.cpp BSDE.cpp
COBJ	= $(CSRC:.c=.o)
CXXOBJ	= $(CXXSRC:.cpp=.o)
OBJ     = $(COBJ) $(CXXOBJ)
CXXOPT  = -O3
LDLIBS  = -lpng
CXXFLAGS= -std=c++0x -fopenmp -Wall
BIN	= BSDE

delault: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(CXXFLAGS) $(LDLIBS) $(CXXOPT)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CXXOPT)

%.o: %.c
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CXXOPT)
clean:
	$(RM) $(OBJ)

distclean: clean
	$(RM) $(BIN) ex*.png

tests: $(BIN)
	@echo Example 1.
	./$(BIN) lena.png 20 ex1_in.png ex1_out.png 1 0.0
	@echo Example 2.
	./$(BIN) lena.png 10 ex2_in.png ex2_out.png 0 0.4
	@echo \n

