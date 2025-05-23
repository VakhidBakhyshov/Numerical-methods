all: a.out

a.out: main.o task_1.o
		g++ main.o task_1.o -o a.out

main.o: main.cpp task_1.hpp
		g++ -c -Wall main.cpp

task_7.o: task_1.cpp task_1.hpp
		g++ -c -Wall task_1.cpp

clean:
		rm -rf *.o a.out