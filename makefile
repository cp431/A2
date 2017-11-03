cancer:	cancer.o
	mpicc cancer.o -o cancer

cancer.o: cancer.c
	mpicc -c -std=c99 cancer.c

clean: 
	rm cancer.o cancer
