sDir = ./src
iDir = ./include
oDir = ./obj

sources  = $(wildcard $(sDir)/*cpp)
includes = $(wildcard $(iDir)/*hpp) $(wildcard $(iDir)/*tpp)
objects  = $(subst $(sDir),$(oDir),$(sources:.cpp=.o))

flags = -I./include

ifeq ($(DOUBLE), 1)
flags += -Ddatafloat=double -Dsize_datafloat=8
else
flags += -Ddatafloat=float -Dsize_datafloat=4
endif

libs = -L. -L./lib -llapack -lblas -lm -ldl -lnetcdf

# OSX
ifeq ($(shell uname -s),Darwin)

libs   += -framework Accelerate -framework CUDA -framework OpenCL -framework OpenGL -framework GLUT
flags += -D 'OS_OSX'
# LINUX
else
libs   += -lcuda -lglut -lGL -lGLU -lGLEW -lOpenCL
flags += -D 'OS_LINUX'
endif


ifeq ($(DEBUG), 1)
# add debugging and force symbol inclusion
flags += -g
LDFLAGS += -rdynamic -g
else
flags += -O3 -DNDEBUG
LDFLAGS += -rdynamic
endif

main: main.cpp $(objects) $(includes)
	$(CXX) $(flags) $(objects) $(libs) main.cpp -o main

$(oDir)/%.o:$(sDir)/%.cpp
	$(CXX) $(flags) -o $@ $(libs) -c $<

clean:
	rm -f $(oDir)/*
	rm -f main
