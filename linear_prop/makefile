SHELL=/bin/bash
DATE = $(shell date +"%d%b%Y")

FC=gfortran
all: qpr.x linden.x
	@echo ""
	@echo "************************************ linden.f **********************************"
	@echo ""
	@echo "[X] Compilation is successful."
	@echo $(shell date)

linden.x:	linden.o
		$(FC) -o $@ linden.o

linden.o:	linden.f
		$(FC) -c linden.f
qpr.x:	qpr.o
		$(FC) -o $@ qpr.o

qpr.o:	qpr.f
		$(FC) -c qpr.f

clean:
		rm -f *.o linden.x
