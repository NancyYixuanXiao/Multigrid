CC= gcc
RM= rm -vf
CFLAGS= -Wall -g -std=c99
OPENMPFLAG= -fopenmp
LIBS= -lrt -lm

PROGFILES= mg

main: serial_mg.o parallel_mg.o helpers.o main.c
	$(CC) $(CFLAGS) $(OPENMPFLAG) serial_mg.o parallel_mg.o helpers.o main.c -o mg $(LIBS)

serial_mg.o: serial_mg.c serial_mg.h
	$(CC) $(CFLAGS) -c serial_mg.c $(LIBS)

parallel_mg.o: parallel_mg.c parallel_mg.h
	$(CC) $(CFLAGS) $(OPENMPFLAG) -c parallel_mg.c $(LIBS)

helpers.o: helpers.c helpers.h
	$(CC) $(CFLAGS) $(OPENMPFLAG) -c helpers.c $(LIBS)

clean:
	$(RM) *.o $(PROGFILES)
