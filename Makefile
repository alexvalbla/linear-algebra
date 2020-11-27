
all: test

matrix.o: matrix.c matrix.h
	gcc -Wall -O3 -c -o matrix.o matrix.c

vector.o: vector.c
	gcc -Wall -O3 -c -o vector.o vector.c

linalg.o: linalg.c
	gcc -Wall -O3 -c -o linalg.o linalg.c

main.o: main.c
	gcc -Wall -c -o main.o main.c

test: matrix.o vector.o linalg.o main.o
	gcc -Wall -o test main.o matrix.o vector.o linalg.o

clean:
	rm *.o test
