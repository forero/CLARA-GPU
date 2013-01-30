
INC	:= -I$(CUDA_HOME)/include -I.
LIB	:= -L$(CUDA_HOME)/lib64 -lcudart

NVCCFLAGS	:= -DCUDA_DEVICE=0 --ptxas-options=-v --use_fast_math

all: 		clara

clara:		main.c init.c scatter.cu 
		nvcc scatter.cu init.c main.c   -o clara \
		$(INC) $(NVCCFLAGS) $(LIB)
