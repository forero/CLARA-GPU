
INC	:= -I$(CUDA_HOME)/include -I.
LIB	:= -L$(CUDA_HOME)/lib64 -lcudart -lm 

NVCCFLAGS	:= -DCUDA_DEVICE=0 --ptxas-options=-v --use_fast_math -m64
EXEC = clara



OBJS = main.o struct.o parse.o io.o init.o transport.o scatter.o
INCLUDES = struct.h parse.h io.h init.h transport.h scatter.h

.SUFFIXES: .o .c 
.c.o:
	nvcc $(NVCCFLAGS) $(INC) $(LIB) -c -lm $<

.SUFFIXES: .o .cu
.cu.o:
	nvcc $(NVCCFLAGS) $(INC) $(LIB) -c -lm $<


all: $(EXEC)

$(EXEC): $(OBJS) $(INCLUDES)
	nvcc $(OBJS) $(INC) $(NVCCFLAGS) $(LIB) -o $(EXEC)

clean:
		rm -f $(EXEC) *.o