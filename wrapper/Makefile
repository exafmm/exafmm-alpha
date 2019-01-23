include ../Makefile

serial: serial.o $(OBJECT)
	$(CXX) $? $(LFLAGS)
	./a.out

lib_serial: serial_wrapper.o $(OBJECT)
	ar -cr libfmm.a serial_wrapper.o $(OBJECT)
	ranlib libfmm.a

test_serial: test_serial.o
	$(FC) $? $(LFLAGS) -L. -lfmm
	./a.out

parallel: parallel.o $(OBJECT)
	$(CXX) $? $(LFLAGS)
	mpirun -np 2 ./a.out

lib_parallel: parallel_wrapper.o $(OBJECT)
	ar -cr libfmm.a parallel_wrapper.o $(OBJECT)
	ranlib libfmm.a

test_parallel: test_parallel.o
	$(FC) $? $(LFLAGS) -L. -lfmm
	mpirun -np 2 ./a.out
