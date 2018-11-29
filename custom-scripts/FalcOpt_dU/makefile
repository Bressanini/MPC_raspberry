


all:
	arm-linux-gcc -c generatedCode/controladorFalcOpt.c 
	arm-linux-gcc -c myBibs/utilidades.c 
	arm-linux-gcc -c mainDeadline.c  
	arm-linux-gcc -o mainDeadline mainDeadline.o controladorFalcOpt.o utilidades.o -lcblas -lm -lpthread
