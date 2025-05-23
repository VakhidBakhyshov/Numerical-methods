all: a.out

a.out: main.o task_3.o
		g++ main.o task_3.o -o a.out

main.o: main.cpp task_3.hpp
		g++ -c -Wall main.cpp

task_3.o: task_3.cpp task_3.hpp
		g++ -c -Wall task_3.cpp

clean:
		rm -rf *.o a.out