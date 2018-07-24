FC=gfortran
all: proposed.x
proposed.x:	proposed.o
		$(FC) -o $@ proposed.o

proposed.o:	proposed.f
		$(FC) -c proposed.f
clean:
		rm -f *.o proposed.x
