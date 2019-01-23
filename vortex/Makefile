include ../Makefile.include

main: main.o $(OBJECT)
	rm -f statistics.dat
	$(CXX) $? $(LFLAGS) -lfftw3
	mpirun -np 2 ./a.out

clean:
	make -C .. clean
