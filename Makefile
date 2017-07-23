CC = g++ -std=gnu++11
CFLAGS = -O3  -fopenmp -Wall -Wextra -Wunused  
#-lboost_thread-mt -O3 -Ofast
LIB = ../../../../lib
DEPS = 
INCLUDES = -I./
GLE_PATH=$(LIB)

OBJ = iview
all:	$(OBJ) mvbin clean

config.o: config.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


iview.o: 	iview.cc 
	$(CC) $(INCLUDES) -c -o $@ $< $(CFLAGS)
iview:	iview.o config.o
	$(CC)  $(CFLAGS) -o $@ $^ -lpng -lz -lgl2ps -lGL -lglut -lX11 -lXext -lGLU -lgle -ljpeg  -lgsl -lgslcblas -lm  -lgomp -lboost_system  -lboost_program_options -lboost_timer  -lboost_chrono


mvbin:
	rm -rf ./bin >& /dev/null; mkdir bin; mv $(OBJ) ./bin
clean:	
	rm -rf *.o *~ $(OBJ) 
src:
	tar -czvf src.tgz *.cc *.h Makefile Doxyfile
ref:
	doxygen Doxyfile; cd doc/latex; make; make
