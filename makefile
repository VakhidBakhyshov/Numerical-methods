all: a.out

a.out: main.o task_2.o
		g++ main.o task_2.o -o a.out

main.o: main.cpp task_2.hpp
		g++ -c -Wall main.cpp

task_2.o: task_2.cpp task_2.hpp
		g++ -c -Wall task_2.cpp

clean:
		rm -rf *.o a.out