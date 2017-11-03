a1:	a1.o
	mpicc a1.o -o a1 -lgmp

a1.o:	a1.c
	mpicc -c -std=c99 a1.c

clean: 
	rm a1.o a1
