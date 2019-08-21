.SUFFIXES: .cxx .cu .o

CXX     = g++ -g -O3 -std=c++11 -fPIC -fopenmp -Iinclude
NVCC    = nvcc -Xcompiler "-fopenmp -O3" -std=c++11 -Iinclude

.cxx.o:
	$(CXX) -c $? -o $@
.cu.o:
	$(NVCC) -c $? -o $@

serial: serial.o gpuLaplace.o
	$(CXX) $? -lcudart
	./a.out

clean:
	@find -name "*.o" -o -name "*.out*" -o -name "*.dat" -o -name "*.a" | xargs rm -rf
