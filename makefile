all: test dotest
test: main.o integrals.o
	g++ -o prog main.o integrals.o --std=c++11 
main.o: main.cpp integrals.h
	g++ -c main.cpp --std=c++11

integrals.o: integrals.cpp integrals.h
	g++ -c integrals.cpp --std=c++11
	
dotest: test
	./prog
