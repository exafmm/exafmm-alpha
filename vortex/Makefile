include ../Makefile.include

main: main.o $(OBJECT)
	rm -f statistics.dat
	$(CXX) $? $(LFLAGS) -lfftw3

clean:
	make -C .. clean
