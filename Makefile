genpotfield: genpotfield.o Makefile
	gcc genpotfield.o -lm -o genpotfield

genpotfield.o:	genpotfield.c Makefile
	gcc -c genpotfield.c
