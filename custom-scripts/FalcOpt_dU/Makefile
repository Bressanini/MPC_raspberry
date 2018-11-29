#include ../lapack-3.7.1/make.inc
include ../lapack-3.7.1/make.inc

CBLASLIB = ../lapack-3.7.1/libcblas.a
BLASLIB  = ../lapack-3.7.1/librefblas.a

FalcOptCod = \
	./generatedCode/casadi_fcn.c \
	./generatedCode/controladorFalcOpt.c

# 	./generatedCode/constraints.c \

all: main

main: main.o $(CBLASLIB) $(BLASLIB) ${FalcOptCod}
	@echo Compilando main
	$(LOADER) $(LOADOPTS) -o $@ $^
	@echo Finalizado

.c.o:
	@echo Compilando main.o
	$(CC) $(CFLAGS) -I../include -c -o $@ $<
	@echo Finalizado

clean: cleanobj cleanexe
cleanobj:
	rm -f *.o
cleanexe:
	rm -f testeMult
	rm -f *.csv
	rm -f sim/X.csv
	rm -f sim/U.csv