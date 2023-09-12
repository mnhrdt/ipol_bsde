all: BSDE

BSDE: BSDE.o io_png.o lib.o 
	g++ -std=c++0x BSDE.o lib.o io_png.o -o BSDE  -fopenmp -lpng -Wall -O3

BSDE.o: BSDE.cpp
	g++ -std=c++0x -c BSDE.cpp -fopenmp -Wall -O3

lib.o: lib.cpp
	g++ -std=c++0x -c lib.cpp -Wall -O3 

io_png.o: io_png.c
	g++ -std=c++0x -c io_png.c -Wall -O3

clean:
	rm -rf *o BSDE
