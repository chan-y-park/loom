CC = gcc
CFLAGS = -Wall -fPIC -std=c99
all: s_wall.so

s_wall.so: s_wall.o solve.o
	$(CC) -Wall -shared -o s_wall.so s_wall.o solve.o

s_wall.o: s_wall.c
	$(CC) $(CFLAGS) -c s_wall.c

solve.o: solve.c
	$(CC) $(CFLAGS) -c solve.c

clean:
	rm *o
