THIS ?= mcmcPaddingtonGAPSV
# THIS ?= particlefilterPaddingtonGAPSV
# note: toolbox2 needs updating before particlefilter can be used

toolboxes=statespacebox.o vslbox.o embox.o gibbsbox.o timerbox.o blaspackbox.o cmcbox.o densitybox.o



# mac: debug
FCdebugseq=ifort -mkl -warn all -WB -check all -check noarg_temp_created -static-intel -Wl,-stack_size,0x20000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/
FCdebug=ifort -mkl -warn all -WB -check all -check noarg_temp_created -static-intel -Wl,-stack_size,0x20000000  -qopenmp -Wl,-rpath,$(MKLROOT)/../compiler/lib/

#  -warn alignments,declarations,general,unused,usage,interfaces

# mac: production
# FC=ifort -O3 -mkl -nocheck -qopenmp -static-intel -xHost -heap-arrays 65532 -Wl,-rpath,$(MKLROOT)/../compiler/lib/
FCprod=ifort -O3 -mkl -nocheck -qopenmp -static-intel -xHost -Wl,-stack_size,0x80000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/

FCprofile=ifort -O3 -mkl -nocheck -qopenmp -static-intel -xHost -qopt-report-file=foo.out -Wl,-stack_size,0x80000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/

# -align array64byte 
# -Wl,-stack_size,0x80000000 -- 2G
# -Wl,-stack_size,0x20000000
# -fast
# -fp-model source

ifeq ($(FCmode),debug)
  FC=$(FCdebug)
else ifeq ($(FCmode),debugseq)
  FC=$(FCdebugseq)
else ifeq ($(FCmode),profile)
  FC=$(FCprofile)
else
  FC=$(FCprod)
endif


toolboxdir=toolbox2/
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

cmcbox.o  : $(toolboxdir)cmcbox.f90 
	$(FC) -c $(toolboxdir)cmcbox.f90 

densitybox.o  : $(toolboxdir)densitybox.f90 embox.o blaspackbox.o vslbox.o statespacebox.o
	$(FC) -c $(toolboxdir)densitybox.f90 

blaspackbox.o  : $(toolboxdir)blaspackbox.f90  embox.o 
	$(FC) -c $(toolboxdir)blaspackbox.f90 

timerbox.o : $(toolboxdir)$(timerbox).f90 embox.o
	$(FC) -c $(toolboxdir)$(timerbox).f90 -o timerbox.o

gibbsbox.o  : $(toolboxdir)gibbsbox.f90 timerbox.o blaspackbox.o statespacebox.o
	$(FC) -c $(toolboxdir)gibbsbox.f90 

vslbox.o  : $(toolboxdir)$(vslbox).f90 
	$(FC) -c $(toolboxdir)$(vslbox).f90 -o vslbox.o

compile	: $(THIS)

run	: $(THIS)
	rm -f *.debug
	time -p ./$(THIS) 


edit :
	aquamacs $(THIS).f90 

clean	:	
	rm -f $(THIS) 
	rm -f *.debug
	rm -f *_genmod.f90
	rm -f *.log
	rm -f *.mod
	rm -f *.o
	rm -f *~

cleanall :
	rm -f *.dat
	$(MAKE) clean	
