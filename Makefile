.SUFFIXES: .cxx .cu .o
.PHONY: docs

CUDA_INSTALL_PATH = /usr/local/cuda

#DEVICE  = cpu
DEVICE  = gpu

CXX     = mpicxx -ggdb3 -O3 -fPIC -fopenmp -ffast-math -funroll-loops -fforce-addr -rdynamic -I../include
NVCC    = nvcc -Xcompiler "-fopenmp -O3" -I../include
LFLAGS  = -D$(DEVICE)
ifeq ($(DEVICE),gpu)
LFLAGS  += -lcudart -lstdc++ -ldl -lm
endif
OBJECT  = ../kernel/$(DEVICE)Laplace.o ../kernel/$(DEVICE)BiotSavart.o\
	../kernel/$(DEVICE)Stretching.o ../kernel/$(DEVICE)Gaussian.o\
	../kernel/$(DEVICE)CoulombVdW.o

.cxx.o:
	$(CXX) -c $? -o $@ $(LFLAGS)

.cu.o:
	$(NVCC) -c $? -o $@ $(LFLAGS)

clean:
	rm -rf `find .. -name "*.o" -o -name "*.out*" -o -name "*.dat" -o -name "*.a"`
