all: s_wall.so

s_wall.so: s_wall.o
	gcc -Wall -shared -o s_wall.so s_wall.o

s_wall.o: s_wall.c
	gcc -c -Wall -fPIC s_wall.c

clean:
	rm *o
