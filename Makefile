all: a.out
a.out: main.o solution.o
	g++ main.o solution.o -o a.out
main.o: main.cpp solution.hpp
	g++ -c  -Wall main.cpp
solution.o: solution.cpp solution.hpp
	g++ -c  -Wall solution.cpp
clean:
	rm -rf *.o a.out
