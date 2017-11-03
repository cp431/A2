a2:	a2.o
	mpicc a2.o -o a2

a2.o: a2.c
	mpicc -c -std=c99 a2.c

clean: 
	rm a2.o a2
