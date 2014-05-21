

all: heat

heat: main.c derivatives.c utilities.c kutta.c change.c
	mpicc -std=c99 -g changeM.c derivatives.c utilities.c Runge.c deriv.c -lm -o heat

