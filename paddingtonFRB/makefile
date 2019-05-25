FORTRANFILE=mcmcPaddingtonGAPSV

FC=ifort -mkl  -warn all -WB -check all -check bounds -check noarg_temp_created

FC=ifort -mkl  -warn all -WB -check all -check bounds -check noarg_temp_created  -openmp # -g -traceback # -diag-enable sc3 -diag-enable warn

FC=ifort -mkl -O3 -nocheck -openmp -inline-level=2 -shared-intel -mcmodel=medium -xSSE4.2 

toolboxdir=toolbox2/
vslbox=INTELvslbox
timerbox=OMPtimerbox

main  : $(FORTRANFILE)

toolboxes : statespacebox.o vslbox.o embox.o gibbsbox.o timerbox.o blaspackbox.o cmcbox.o densitybox.o

$(FORTRANFILE)  : $(FORTRANFILE).f90 statespacebox.o vslbox.o embox.o gibbsbox.o timerbox.o blaspackbox.o cmcbox.o densitybox.o
	$(FC) $(FORTRANFILE).f90 -o $(FORTRANFILE) embox.o vslbox.o gibbsbox.o statespacebox.o timerbox.o blaspackbox.o cmcbox.o densitybox.o

statespacebox.o :: $(toolboxdir)statespacebox.f90 blaspackbox.o embox.o vslbox.o
	$(FC) -c $(toolboxdir)statespacebox.f90

embox.o  : $(toolboxdir)embox.f90 vslbox.o
	$(FC) -c $(toolboxdir)embox.f90 

cmcbox.o  : $(toolboxdir)cmcbox.f90 
	$(FC) -c $(toolboxdir)cmcbox.f90 

densitybox.o  : $(toolboxdir)densitybox.f90 embox.o blaspackbox.o vslbox.o
	$(FC) -c $(toolboxdir)densitybox.f90 

blaspackbox.o  : $(toolboxdir)blaspackbox.f90  embox.o 
	$(FC) -c $(toolboxdir)blaspackbox.f90 

timerbox.o : $(toolboxdir)$(timerbox).f90 embox.o
	$(FC) -c $(toolboxdir)$(timerbox).f90 -o timerbox.o

gibbsbox.o  : $(toolboxdir)gibbsbox.f90 timerbox.o blaspackbox.o 
	$(FC) -c $(toolboxdir)gibbsbox.f90 

vslbox.o  : $(toolboxdir)$(vslbox).f90 
	$(FC) -c $(toolboxdir)$(vslbox).f90 -o vslbox.o

run	: $(FORTRANFILE)
	rm -f *.debug
	time -p ./$(FORTRANFILE) 

clean	:	
	rm -f $(FORTRANFILE) 
	rm -f *.dat
	rm -f *.debug
	rm -f *_genmod.f90
	rm -f *.log
	rm -f *.mod
	rm -f *.o
	rm -f *.o*
	rm -f *.e*
	rm -f *~

