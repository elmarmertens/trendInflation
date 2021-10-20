THIS ?= mcmcPaddingtonGAPSV
datalabel ?= INFTRM
# THIS ?= particlefilterPaddingtonGAPSV
# note: toolbox2 needs updating before particlefilter can be used

toolboxes=statespacebox.o vslbox.o embox.o gibbsbox.o timerbox.o blaspackbox.o densitybox.o

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
  # mac
  FCfulldebugseq=ifort -qmkl -warn all -warn noexternals -WB -check all -check bounds -g -check noarg_temp_created -static-intel -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -Wl,-stack_size,0x20000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  FCdebugseq=ifort -qmkl -warn all -warn noexternals -WB -check all -check noarg_temp_created -static-intel -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -Wl,-stack_size,0x20000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  FCfulldebug=ifort -qmkl -warn all -WB -warn noexternals -check all -check bounds -traceback -g -check noarg_temp_created -static-intel -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -Wl,-stack_size,0x20000000  -qopenmp -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  FCdebug=ifort -qmkl -warn all -WB -warn noexternals -check all -check noarg_temp_created -static-intel -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -Wl,-stack_size,0x20000000  -qopenmp -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  FCprod=ifort -O3 -qmkl -nocheck -qopenmp -static-intel -xHost -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -Wl,-stack_size,0x80000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/
else
  # linux
  FCdebugseq=ifort -qmkl -warn all -WB -check all -check noarg_temp_created  -shared-intel 
  FCdebug=ifort -qmkl -warn all -WB -check all -check noarg_temp_created  -shared-intel -qopenmp 
  FCfulldebug=ifort -qmkl -warn all -WB -check all -check noarg_temp_created  -shared-intel -qopenmp 
  FCprod=ifort -O3 -qmkl -nocheck -qopenmp  -shared-intel  -xHost 
endif

ifeq ($(FCmode),debug)
  FC=$(FCdebug)
else ifeq ($(FCmode),debugseq)
  FC=$(FCdebugseq)
else ifeq ($(FCmode),fulldebug)
  FC=$(FCfulldebug)
else ifeq ($(FCmode),fulldebugseq)
  FC=$(FCfulldebugseq)
else
  FC=$(FCprod)
endif


toolboxdir=fortranbox/
vslbox=INTELvslbox
timerbox=OMPtimerbox

main  : $(THIS)
toolboxes : $(toolboxes)

$(THIS)  : $(THIS).f90 $(toolboxes) 
	$(FC) $(THIS).f90 -o $(THIS) $(toolboxes) 



statespacebox.o :: $(toolboxdir)statespacebox.f90 blaspackbox.o embox.o vslbox.o
	$(FC) -c $(toolboxdir)statespacebox.f90

embox.o  : $(toolboxdir)embox.f90 vslbox.o
	$(FC) -c $(toolboxdir)embox.f90 

densitybox.o  : $(toolboxdir)densitybox.f90 embox.o blaspackbox.o vslbox.o statespacebox.o
	$(FC) -c $(toolboxdir)densitybox.f90 

blaspackbox.o  : $(toolboxdir)blaspackbox.f90  embox.o 
	$(FC) -c $(toolboxdir)blaspackbox.f90 

timerbox.o : $(toolboxdir)$(timerbox).f90 embox.o
	$(FC) -c $(toolboxdir)$(timerbox).f90 -o timerbox.o

gibbsbox.o  : $(toolboxdir)gibbsbox.f90 timerbox.o blaspackbox.o statespacebox.o densitybox.o
	$(FC) -c $(toolboxdir)gibbsbox.f90 

vslbox.o  : $(toolboxdir)$(vslbox).f90 
	$(FC) -c $(toolboxdir)$(vslbox).f90 -o vslbox.o

compile	: $(THIS) 

run	: $(THIS)
	rm -f *.debug
	time -p ./$(THIS) $(datalabel)

edit : 
	aquamacs $(THIS).f90 

clean	:	
	rm -f $(THIS) 
	rm -f *.debug
	rm -rf *.dSYM
	rm -f *_genmod.f90
	rm -f *.log
	rm -f *.mod
	rm -f *.o
	rm -f *~

cleanall :
	rm -f *.dat
	rm -f *.jpg
	rm -f *.csv
	$(MAKE) clean	
