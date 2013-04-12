CC=cc -D SKIP_CLOCK
CFLAGS=-O3 -fopenmp -fnested-functions
LIBS=-lm -lgsl -lgslcblas -L/opt/local/lib -I/opt/local/include
RM=/bin/rm

EXECS = ko
OBJS = ko.o fileop.o misc.o physics.o finite.o zsolve_quartic.o problem.o metric.o relele.o rad.o rad.mf.o u2p.o frames.o p2u.o

all: $(EXECS)

$(EXECS): $(OBJS) Makefile ko.h problem.h mnemonics.h Makefile.gcc Makefile.icc
	$(CC) $(CFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

$(OBJS): Makefile ko.h problem.h mnemonics.h Makefile.gcc Makefile.icc

.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c
clean:
	$(RM) -f $(EXECS) $(OBJS) *~ *.o *.oo
