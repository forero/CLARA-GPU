
INC	:= -I$(CUDA_HOME)/include -I.
LIB	:= -L$(CUDA_HOME)/lib64 -lcudart

NVCCFLAGS	:= -DCUDA_DEVICE=0 --ptxas-options=-v --use_fast_math

all: 		clara

clara:		main.c io.c init.c transport.c scatter.cu 
		nvcc scatter.cu transport.c init.c io.c main.c   -o clara \
		$(INC) $(NVCCFLAGS) $(LIB)
