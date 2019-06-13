include ../Makefile.include

serial: serial.o $(OBJECT)
	$(CXX) $? $(LFLAGS)

lib_serial: serial_wrapper.o $(OBJECT)
	ar -cr libfmm.a serial_wrapper.o $(OBJECT)
	ranlib libfmm.a

test_serial: test_serial.o
	$(FC) $? $(LFLAGS) -L. -lfmm

parallel: parallel.o $(OBJECT)
	$(CXX) $? $(LFLAGS)

lib_parallel: parallel_wrapper.o $(OBJECT)
	ar -cr libfmm.a parallel_wrapper.o $(OBJECT)
	ranlib libfmm.a

test_parallel: test_parallel.o
	$(FC) $? $(LFLAGS) -L. -lfmm

clean:
	make -C .. clean
