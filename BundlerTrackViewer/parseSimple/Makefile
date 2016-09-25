CC        := g++
LD        := g++

MODULES   := .
SRC_DIR   := $(addprefix ./, $(MODULES))
BUILD_DIR := $(addprefix ./, $(MODULES))
SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst %.cpp,%.o,$(SRC))

YOURINCLUDEPATH = -I/Users/ceylan/MVSReconstruction/bundlerSource/lib/matrix -I/Users/ceylan/OpenMesh/src

# Your LDLIBRARY path (e.g., -L/usr/lib)
YOURLDLIBPATH = -L/usr/local/lib -L/Users/ceylan/MVSReconstruction/bundlersource/lib/matrix -L/Users/ceylan/MVSReconstruction/bundlerSource/lib

CXXFLAGS = ${YOURINCLUDEPATH} ${INCLUDES}

LDFLAGS = ${YOURLDLIBPATH} -lmatrix -ljpeg -lpng -lz -lminpack -framework Accelerate

define make-goal
$1/%.o: %.cpp
	$(CC) -c $(CXXFLAGS) -c $$< -o $$@
endef

all: parse

parse: Image.o parseBundlerOutput.o  $(OBJ)
	$(LD) ${LDFLAGS} -o $@ $^ ${LDFLAGS}

Image.o: Image.cpp 
	$(CC) -c $(CXXFLAGS) $<
parseBundlerOutput.o: parseBundlerOutput.cpp
	$(CC) -c $(CXXFLAGS) $<
clean:
	rm -rf $(OBJ); rm main.o     