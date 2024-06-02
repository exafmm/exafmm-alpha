include ../Makefile.include

serial: serial.o $(OBJECT)
	$(CXX) $? $(LFLAGS)
	./a.out

lib_serial: serial_wrapper.o $(OBJECT)
	rm -f libfmm.a
	ar -cr libfmm.a serial_wrapper.o $(OBJECT)
	ranlib libfmm.a

test_serial: test_serial.o
	make lib_serial
	$(FC) $? -L. -lfmm $(LFLAGS)
	./a.out

lib_serial_ij: serial_wrapper_ij.o $(OBJECT)
	rm -f libfmm.a
	ar -cr libfmm.a serial_wrapper_ij.o $(OBJECT)
	ranlib libfmm.a

test_serial_ij: test_serial_ij.o
	make lib_serial_ij
	$(FC) $? -L. -lfmm $(LFLAGS)
	./a.out

parallel: parallel.o $(OBJECT)
	$(CXX) $? $(LFLAGS)
	mpirun -np 2 ./a.out

lib_parallel: parallel_wrapper.o $(OBJECT)
	rm -f libfmm.a
	ar -cr libfmm.a parallel_wrapper.o $(OBJECT)
	ranlib libfmm.a

test_parallel: test_parallel.o
	make lib_parallel
	$(FC) $? -L. -lfmm $(LFLAGS)
	mpirun -np 2 ./a.out

lib_parallel_ij: parallel_wrapper_ij.o $(OBJECT)
	rm -f libfmm.a
	ar -cr libfmm.a parallel_wrapper_ij.o $(OBJECT)
	ranlib libfmm.a

test_parallel_ij: test_parallel_ij.o
	make lib_parallel_ij
	$(FC) $? -L. -lfmm $(LFLAGS)
	mpirun -np 2 ./a.out

clean:
	make -C .. clean
