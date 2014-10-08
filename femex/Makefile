CXX = g++

CXX_FLAGS = -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
			-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread\
			-DMX_COMPAT_32 -O -DNDEBUG

CXX_INCLUDE = -I./include/mexplus \
			  -I./include/pprint \
			  -I/usr/local/MATLAB/MATLAB_Production_Server/R2013a/extern/include \
			  -I/usr/local/MATLAB/MATLAB_Production_Server/R2013a/simulink/include


MATLAB_LINKS = -O -pthread -shared\
			   -Wl,--version-script,/usr/local/MATLAB/MATLAB_Production_Server/R2013a/extern/lib/glnxa64/mexFunction.map \
			   -Wl,--no-undefined 
			   
CXX_LIBS = -Wl,-rpath-link,/usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/glnxa64 \
		   -L/usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/glnxa64 -lmx -lmex -lmat -lm
		   
		   
		   
SRC = src/
# mesh 
MESH = $(SRC)Mesh/private/
TRIANGLELIB = $(MESH)triangle/
MESH_SRCS = $(wildcard $(MESH)*.cc)
MESH_OBJS = $(patsubst $(MESH)%.cc, %.o, $(MESH_SRCS))
MESH_BINS = $(patsubst $(MESH)%.cc, $(MESH)%_.mexa64, $(MESH_SRCS))

%.o: $(MESH)%.cc
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(MESH)%_.mexa64: %.o $(TRIANGLELIB)triangle.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(TRIANGLELIB)triangle.o $(CXX_LIBS)
	
	
# Triangle build

CSWITCHES = -O -DLINUX -fPIC -shared -I/usr/X11R6/include -L/usr/X11R6/lib
TRILIBDEFS = -DTRILIBRARY

$(TRIANGLELIB)triangle.o:$(TRIANGLELIB)triangle.c $(TRIANGLELIB)triangle.h
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -c -o $(TRIANGLELIB)triangle.o \
	$(TRIANGLELIB)triangle.c 



# The action starts here.

# all: $(BIN)MeshGen.mexa64 $(SRC)triangle.o $(BIN)MatrixAssem.mexa64 $(BIN)BCAssem.mexa64 clean

all: $(MESH_BINS)

clean:
	rm -rf $(MESH)*_.mexa64 