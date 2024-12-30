ARGS = 20

.PHONY: all clean run graphs

all: a.out

a.out: main.o task_5.o
		g++ main.o task_5.o  -o a.out

main.o: main.cpp task_5.hpp
		g++ -c -Wall main.cpp

task_5.o: task_5.cpp task_5.hpp
		g++ -c -Wall task_5.cpp

run: a.out
		./a.out $(ARGS) 0 out.txt

clean:
		rm -rf *.o a.out
