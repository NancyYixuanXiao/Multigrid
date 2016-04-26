CC= gcc
RM= rm -vf
CFLAGS= -Wall -g
OPENMPFLAG= -fopenmp
LIBS= -lm -lrt

PROGFILES= mg

main: serial_mg.o parallel_mg.o main.c
	$(CC) $(CFLAGS) $(LIBS) $(OPENMPFLAG) serial_mg.o parallel_mg.o main.c -o mg

serial_mg.o: serial_mg.c serial_mg.h
	$(CC) $(CFLAGS) $(LIBS) -c serial_mg.c

parallel_mg.o: parallel_mg.c parallel_mg.h
	$(CC) $(CFLAGS) $(LIBS) $(OPENMPFLAG) -c parallel_mg.c

clean:
	$(RM) *.o $(PROGFILES)
