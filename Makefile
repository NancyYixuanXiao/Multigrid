CC= gcc
RM= rm -vf
CFLAGS= -Wall -g
OPENMPFLAG= -fopenmp
LIBS= -lm -lrt

PROGFILES= mg

main: mg.c
	$(CC) $(CFLAGS) $(LIBS) $(OPENMPFLAG) mg.c -o mg

clean:
	$(RM) *.o $(PROGFILES)
