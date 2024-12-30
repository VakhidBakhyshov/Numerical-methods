ARGS = 20

all: a.out

a.out: main.o task_1.o pcalculate.o
	g++ main.o task_1.o pcalculate.o -o a.out

main.o: main.cpp task_1.hpp
	g++ -c -Wall main.cpp

task_1.o: task_1.cpp task_1.hpp
	g++ -c -Wall task_1.cpp

pcalculate.o: pcalculate.cpp
	g++ -c -Wall pcalculate.cpp

run: a.out
	./a.out $(ARGS) 0 out.txt

clean:
	rm -rf *.o a.out
