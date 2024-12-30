ARGS = 20

.PHONY: all clean run graphs

all: a.out

a.out: main.o task_7.o
	g++ main.o task_7.o  -o a.out

main.o: main.cpp task_7.hpp
	g++ -c -Wall main.cpp

task_7.o: task_7.cpp task_7.hpp
	g++ -c -Wall task_7.cpp

clean:
	rm -rf *.o a.out
