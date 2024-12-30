ARGS = 20
.PHONY: all clean run

all: a.out

a.out: main.o task_4.o system.o
		g++ $^ -o $@

main.o: main.cpp task_4.hpp system.hpp
		g++ -c -Wall main.cpp

task_4.o: task_4.cpp task_4.hpp system.hpp
		g++ -c -Wall task_4.cpp

system.o: system.cpp system.hpp
		g++ -c -Wall system.cpp

run: a.out
		./a.out $(ARGS) 0 out.txt

clean:
		rm -rf *.o a.out