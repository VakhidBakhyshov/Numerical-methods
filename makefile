ARGS = 20

.PHONY: all clean run

all: a.out

a.out: main.o task_3.o system.o
	g++ $^ -o $@

main.o: main.cpp task_3.hpp system.hpp
	g++ -c -Wall main.cpp

task_3.o: task_3.cpp task_3.hpp system.hpp
	g++ -c -Wall task_3.cpp

system.o: system.cpp system.hpp
	g++ -c -Wall system.cpp

run: a.out
	./a.out $(ARGS) 0 out.txt

clean:
	rm -rf *.o a.out
