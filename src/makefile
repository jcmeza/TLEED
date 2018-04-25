#
#	Create TLEED/KLEED Standalone Executables
#
###############################################################
# gcc compilers
FC     = gfortran 
FFLAGS  = -O2 -Wall -Wsurprising -Waliasing  -fbacktrace -fdefault-real-8
#FFLAGS  = -O2 -Wall -Wsurprising -Waliasing  -fbacktrace -Wuninitialized -fdefault-real-8

# Intel compiler
#FC     = ifort 
#FFLAGS	= -O -check bounds -check stack -traceback -autodouble  	
#FFLAGS	= -O2 -autodouble  	

FLINKER = $(FC)
###############################################################

SOURCES = evaltleed.f evalkleed.f tleed1.f tleed2.f tleedlib.f kleed.f tleed1sub.f tleed2sub.f kleedsub.f kleedlib.f
EXECS   = tleed1 tleed2 kleed runtleed.exe runkleed.exe kleedfcn.exe

%.o:	%.f	
	$(FC) $(FFLAGS) -c $< -o $@

all:	$(EXECS)

tleed1:	tleed1.o tleedlib.o
	$(FC) $(FFLAGS) -o tleed1.exe tleed1.o tleedlib.o

tleed2:	tleed2.o tleedlib.o
	$(FC) $(FFLAGS) -o tleed2.exe tleed2.o tleedlib.o

kleed:	kleed.o kleedlib.o
	$(FC) $(FFLAGS) -o kleed.exe kleed.o kleedlib.o

runtleed.exe: runtleed.o evaltleed.o tleed1sub.o tleed2sub.o tleedlib.o
	$(FC) $(FFLAGS) -o runtleed.exe runtleed.o evaltleed.o tleed1sub.o tleed2sub.o tleedlib.o

runkleed.exe: runkleed.o evalkleed.o kleedsub.o kleedlib.o
	$(FC) $(FFLAGS) -o runkleed.exe runkleed.o evalkleed.o kleedsub.o kleedlib.o
#
# Matlab/KLEED interface routine
kleedfcn.exe: kleedfcn.o evalkleed.o kleedsub.o kleedlib.o
	$(FC) $(FFLAGS) -o kleedfcn.exe kleedfcn.o evalkleed.o kleedsub.o kleedlib.o

nomad:
	make -f makenomadtleed

clean:
	rm -f *.o *.exe *~ *.f90 *.mod
