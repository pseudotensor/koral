CC=gcc 
CFLAGS=-g -O3 -fopenmp -L/usr/local/lib -I/usr/local/include 
LIBS=-lm -lgsl -lgslcblas -lsilo
RM=/bin/rm

OBJS = u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o zsolve_quartic.o problem.o metric.o relele.o rad.o rad.mf.o u2p.o frames.o p2u.o

all: ko ana

ko: ko.o $(OBJS) Makefile ko.h problem.h mnemonics.h Makefile.gcc Makefile.icc
	$(CC) $(CFLAGS) -o ko ko.o $(OBJS) $(LIBS)

ana: ana.o $(OBJS)  Makefile ko.h problem.h mnemonics.h Makefile.gcc Makefile.icc
	$(CC) $(CFLAGS) -o ana ana.o $(OBJS) $(LIBS)

clean:
	$(RM) -f ko ana $(OBJS) *~ *.o *.oo
