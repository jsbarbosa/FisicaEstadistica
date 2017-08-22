all : *.dat plots.py
	python plots.py

*.dat : a.out
	./a.out

a.out : randomWalk.c
	gcc randomWalk.c -lm -O2 -Wall
