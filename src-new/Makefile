CUDA_HOME=/usr/local/cuda/
INC	:= -I$(CUDA_HOME)/include -I.
LIB	:= -L$(CUDA_HOME)/lib64 -lcudart -lm 

NVCCFLAGS	:= -DDEBUG -DCUDA_DEVICE=0 --ptxas-options=-v --use_fast_math -m64
EXEC = clara

OBJS = main.o struct.o parse.o io.o myrand.c ray.c propagation.o propagation_kernel.o mtwist.o
INCLUDES = struct.h parse.h io.h mtwist.h 

.SUFFIXES: .o .c 
.c.o:
	nvcc $(NVCCFLAGS) $(INC) $(LIB) -c -lm $<
#	cc $(INC) $(LIB) -c -lm $<

.SUFFIXES: .o .cu
.cu.o:
	nvcc $(NVCCFLAGS) $(INC) $(LIB) -c -lm $<


all: $(EXEC)

$(EXEC): $(OBJS) $(INCLUDES)
	nvcc $(OBJS) $(INC) $(NVCCFLAGS) $(LIB) -o $(EXEC)

clean:
		rm -f $(EXEC) *.o