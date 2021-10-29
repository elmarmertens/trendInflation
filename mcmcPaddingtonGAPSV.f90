PROGRAM main

  ! command line arguments: datalabel, Ny, Nsim, burnin, Nstreams, p, sqrtVf0_firstlag, sqrtVf0_general

  USE embox, only : hrulefill, savemat, savevec, storeEstimates, storeEstimatesOMP, loft, timestampstr, es30d16, int2str
  USE blaspack, only : eye, vech
  USE gibbsbox, only : igammaDraw, iwishDraw, GelmanTest1, simpriormaxroot
  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  INTEGER :: Ny, p, Nbar, Ngap, Nstates, Nx, Nsv, NhSigma, Nf, Nshockslopes

  LOGICAL, PARAMETER :: doTimestamp = .false., doDiffuse = .false.

  double precision, parameter :: shockslopesSTD = 1.0d-1

  double precision :: sqrtVf0_firstlag, sqrtVf0_general

  INTEGER :: Nsim, Burnin, Nstreams
  INTEGER :: T,j,k,i,status,Ndraws

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWsvol, DRAWstates
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWf, DRAWshockslopes, DRAWmaxlambda, DRAWhSigma, DRAWhvarbar, DRAWhgap0

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E0shockslopes, Eh0, minSV, maxSV, Ef0
  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: hSigmaT, sqrtVh0
  INTEGER :: hSigmaDof, hvarbarDof
  DOUBLE PRECISION :: hvarbarT, Ehbar0, Vhbar0

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SigmaY, sqrtVf0, sqrtV0shockslopes

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: GELMANstates, GELMANsvol, theta
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GELMANf, GELMANshockslopes, GELMANhSigma, GELMANhgap0, GELMANhvarbar

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ydata
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN


  TYPE(progresstimer) :: timer
  CHARACTER (LEN=200) :: filename, datafile, nandatafile, filext, datalabel

  ! VSL Random Stuff
  type (vsl_stream_state) :: VSLdefaultstream, VSLstream
  integer :: seed
  integer :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OPEN MP
  INTEGER :: NTHREADS, TID

  ! ----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  ! ----------------------------------------------------------------------------------
  ! runtime parameters :start:
  ! first: set default values

  Nsim   = 2 * (10 ** 4)
  Burnin = 10 ** 6

  Nsim    = 10 ** 4
  Burnin  = 10 ** 4

  Nsim    = 10 ** 3
  Burnin  = 10 ** 3

  datalabel = 'INFTRM'

  ! original settings
  sqrtVf0_firstlag = 0.10d0
  sqrtVf0_general  = 0.01d0


  p = 12

  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  call hrulefill
  print *, 'MCMC estimation of PaddingtonGAPSV'
  call hrulefill


  ! INIT OMP
  NTHREADS = 1
  !$OMP PARALLEL SHARED(NTHREADS)
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  Nstreams = max(NTHREADS, 4)
  print *, "Number of Threads:", NTHREADS
  print *, "Number of Streams:", Nstreams

  ! VSL
  call hrulefill
  write (*,*) 'Allocating VSLstreams ...'

  seed    = 0

  errcode = vslnewstream(VSLdefaultstream, vsl_brng_mt2203, seed)  
  WRITE(*,*) "... default VSLstream ", VSLdefaultstream
  T = 0
  call getarguments(datalabel,p,T)
  call getsettings(datalabel,Ny)

  call hrulefill
  print *, 'datalabel: ', datalabel
  if (doDiffuse) then
     print *, 'with diffuse prior'
  else
     print *, 'with tight prior'
  end if
  call hrulefill

  datafile   = trim(datalabel) // '.yData.txt'

  if (T == 0) then
     ! read data
     T = loft(datafile) 
     if (T < 10) THEN
        print *, 'Less than 10 observations in input file!', datafile
        stop 1
     end IF
  end if

  ! runtime parameters :end: 

  filext = '.' // trim(datalabel) // '.T' // trim(int2str(T)) // '.gapSV.dat'
  if (doDiffuse) filext = '.diffuse' //  filext
  IF (doTimeStamp) filext = '.' // timestampstr() // filext
  filext = '.' // 'notrendslopes' // filext


  Nbar = Ny
  Ngap = Ny 
  Nstates = Nbar + Ngap
  Nx = Nbar + p * Ngap
  Nsv = 1 + Ny
  NhSigma = (Ngap + 1) * Ngap / 2
  Nshockslopes = (Ngap - 1) * Ngap / 2
  Nf  = Ngap * Ngap * p


  datafile    = trim(datalabel) // '.yData.txt'
  nandatafile = trim(datalabel) // '.yNaN.txt'

  ALLOCATE (yNaN(Ny,T), ydata(Ny,T), SigmaY(Ny,Ny), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Ydata)'
  END IF

  CALL readdata(ydata,datafile,Ny,T)
  CALL readnandata(yNaN,nandatafile,Ny,T)


  ! print *, 'START ALLOCATION'
  ! print *, 'Ny', Ny
  ! print *, 'T', T
  ! print *, 'p', p
  ! print *, 'Nf', Nf

  ! Model parameters and priors
  ALLOCATE (E0shockslopes(Nshockslopes), sqrtV0shockslopes(Nshockslopes,Nshockslopes), hSigmaT(Nsv-1,Nsv-1), Eh0(Nsv-1), sqrtVh0(Nsv-1,Nsv-1), minSV(Nsv), maxSV(Nsv), Ef0(Nf), sqrtVf0(Nf,Nf))

  print *, 'DONE ALLOCATION'


  ! SV
  minSV = 1.0d-3
  maxSV = 40.0d0 

  ! TREND SV
  Vhbar0     = 1.0d0
  Ehbar0     = log(1 / 12.0d0) - Vhbar0 * 0.25d0;
  hvarbarDof = 2 + 24
  hvarbarT   = (0.2d0 ** 2) / 3.0d0 * dble(hvarbarDof - 2)


  ! Gap SV
  Eh0        = log(0.5d0 ** 2) ! 0.0d0
  call eye(sqrtVh0, 1.0d0)

  Eh0        = log(0.25d0 ** 2) ! 0.0d0
  call eye(sqrtVh0, 2.0d0)

  hSigmaDof = Ngap + 1 + 24
  call eye(hSigmaT, ((.2d0 / sqrt(3.0d0)) ** 2) * dble(hSigmaDof - Ngap - 1))

  ! VAR coefficients
  Ef0  = 0.0d0
  call eye(sqrtVf0, sqrtVf0_general)
  ! widen the prior for first-own-lag coefficients
  DO j=1,NGap
     k = (j-1)*(Ngap * p) + j
     sqrtVf0(k,k) = sqrtVf0_firstlag
  END DO

  ! shockslopes
  E0shockslopes = 0.0d0
  call eye(sqrtV0shockslopes, shockslopesSTD)

  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data= ' // datalabel
  print *, 'T= ', T
  print *, 'Ny= ', Ny
  print *, 'Nsim= ', Nsim
  print *, 'Burnin= ', Burnin
  print *, 'Nstreams= ', Nstreams
  print *, 'p= ', p
  print *, 'sqrtVf0_firstlag= ', sqrtVf0_firstlag
  print *, 'sqrtVf0_general = ', sqrtVf0_general
  print *, 'sqrtV0shockslopes= ', sqrtV0shockslopes(1,1)
  CALL HRULEFILL

  ! allocate memory for draws
  ALLOCATE (DRAWmaxlambda(1,Nsim,Nstreams), DRAWhSigma(NhSigma, Nsim, Nstreams), DRAWhgap0(Ngap, Nsim, Nstreams), DRAWhvarbar(1,Nsim, Nstreams), DRAWf(Nf,Nsim,Nstreams), DRAWshockslopes(Nshockslopes,Nsim,Nstreams), DRAWstates(Nstates,0:T,Nsim,Nstreams), DRAWsvol(Nsv,0:T,Nsim,Nstreams), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (draws 1)'
  END IF

  WRITE(*,*) 'STARTING SWEEPS'

  !$OMP PARALLEL DO SHARED(DRAWstates,DRAWsvol,DRAWf,DRAWmaxlambda,Nsim,Burnin) PRIVATE(TID,timer,VSLstream)
  DO j=1,Nstreams

     TID = 0
     !$ TID = OMP_GET_THREAD_NUM()

     errcode = vslnewstream(VSLstream, vsl_brng_mt2203 + 1 + j, 0)  
     if (errcode /= 0) then
        print *,'VSL new stream failed'
        stop 1
     end if
     WRITE(*,'(a16, i2, a5, i2, a8, i20, i20)') 'LAUNCHING SWEEP ', j, ' TID ', TID, ' Stream ', VSLstream%descriptor1, VSLstream%descriptor2





     WRITE(*,*) 'LAUNCHING SWEEP ', j
     CALL initprogressbar(timer, 15.0d0, j)
     CALL thissampler(doDiffuse,T,p,ydata,yNaN,Ny,DRAWstates(:,:,:,j),Nstates,Nx,DRAWsvol(:,:,:,j),Nsv,Ehbar0, Vhbar0, Eh0,sqrtVh0,minSV,maxSV, DRAWhvarbar(:,:,j), hvarbarT, hvarbarDof, DRAWhgap0(:,:,j), DRAWhSigma(:,:,j), NhSigma, hSigmaT, hSigmaDof, DRAWshockslopes(:,:,j), Nshockslopes, E0shockslopes, sqrtV0shockslopes, DRAWf(:,:,j), Nf, Ef0, sqrtVf0, DRAWmaxlambda(:,:,j), Nsim,Burnin, VSLstream,timer)

     WRITE(*,*) 'STREAM', j, 'IS DONE.', ' (TID: ', TID, ')'
     errcode = vsldeletestream(VSLstream)     

  END DO
  !$OMP END PARALLEL DO 


  CALL HRULEFILL
  WRITE (*,*) 'ALL STREAMS ARE DONE!'
  CALL HRULEFILL

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a20)') 'TIME: ', timestampstr()
  WRITE(4,'(a20,a20)') 'Data: ', datalabel
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a20,I20)') 'Sims:', Nsim
  WRITE(4,'(a20,I20)') 'Burnin:', Burnin
  WRITE(4,'(a20,I20)') 'Streams:', Nstreams
  WRITE(4,'(a20,I20)')   'p               = ', p
  if (doDiffuse) then
     WRITE(4,'(a20)') 'with diffuse VAR prior'
  else
     WRITE(4,'(a20,f10.3)') 'sqrtVf0_firstlag= ', sqrtVf0_firstlag
     WRITE(4,'(a20,f10.3)') 'sqrtVf0_general = ', sqrtVf0_general
     WRITE(4,'(a20,f10.3)') 'sqrtV0shockslopes:', sqrtV0shockslopes(1,1)
  end if
  WRITE(4,'(a20,f10.3)') 'trend min SV:', minSV(1)
  WRITE(4,'(a20,f10.3)') 'gap   min SV:', minSV(2)
  CLOSE(UNIT=4)
  CALL HRULEFILL


  ! ----------------------------------------------------------------------------
  ! STORE
  ! ----------------------------------------------------------------------------
  CALL HRULEFILL
  WRITE (*,*) 'STARTING W/STORAGE'
  CALL HRULEFILL

  ! STORE ESTIMATES
  ! Note: manual reshape avoids segmentation faults
  Ndraws = Nsim * Nstreams

  filename = 'YDATA' // filext
  call savemat(ydata, filename)

  ALLOCATE (theta(T,Ndraws))

  DO i=1,Nstates
     filename  = 'STATES' // trim(int2str(i)) // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(i,1:T,k,j)
     CALL storeEstimates(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED STATES', i
  END DO

  DO i=1,Ny
     filename  = 'TAU' // trim(int2str(i)) // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(i,1:T,k,j)
     CALL storeEstimates(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED TAU', i
  END DO

  filename  = 'SIGTREND' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWsvol(1,1:T,k,j)
  CALL storeEstimates(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED SIGTREND'

  DO i=1,Ny
     filename  = 'SIGGAP' // trim(int2str(i)) // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWsvol(1+i,1:T,k,j)
     CALL storeEstimates(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED SIGGAP', i
  END DO

  DEALLOCATE(theta)

  ALLOCATE (theta(Nbar,Ndraws))
  filename  = 'INITIALTAU' // filext
  FORALL (i=1:Nbar,k=1:Nsim,j=1:Nstreams) theta(i,(j-1) * Nsim + k) = DRAWstates(i,1,k,j)
  ! call savemat(transpose(theta), filename)
  CALL storeEstimates(theta,Nbar,Ndraws,filename)
  WRITE (*,*) 'STORED INITIAL TAU'
  DEALLOCATE(theta)



  ALLOCATE(theta(NhSigma,Ndraws))
  filename  = 'HSIGMA.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhSigma(:,k,j)
  CALL savemat(transpose(theta), filename)
  filename  = 'HSIGMA' // filext
  CALL storeEstimates(theta,NhSigma,Ndraws,filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HSIGMA'

  ALLOCATE(theta(Ngap,Ndraws))
  filename  = 'HGAP0' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhgap0(:,k,j)
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HGAP0'

  ALLOCATE(theta(1,Ndraws))
  filename  = 'HVARBAR' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhvarbar(:,k,j)
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HVARBAR'

  ALLOCATE(theta(Nf,Ndraws))
  filename  = 'F' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWf(:,k,j)
  CALL storeEstimates(theta,Nf,Ndraws,filename)
  ! CALL savemat(theta, filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED F'

  ALLOCATE(theta(Nshockslopes,Ndraws))
  filename  = 'SHOCKSLOPES.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWshockslopes(:,k,j)
  CALL savemat(transpose(theta), filename)
  filename  = 'SHOCKSLOPES' // filext
  CALL storeEstimates(theta,Nshockslopes,Ndraws,filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED SHOCKSLOPES'

  ALLOCATE(theta(1,Ndraws))
  filename  = 'MAXLAMBDA' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWmaxlambda(:,k,j)
  ! CALL storeEstimates(theta,1,Ndraws,filename)
  CALL savevec(theta(1,:), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED MAXLAMBDA'


  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  ! GELMAN
  ! ----------------------------------------------------------------------------

  ALLOCATE (GELMANstates(Nstates,T), GELMANsvol(Nsv,T), GELMANhgap0(Ngap), GELMANhSigma(NhSigma), GELMANhvarbar(1), GELMANf(Nf), GELMANshockslopes(Nshockslopes), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Gelman statistics)'
  END IF

  CALL HRULEFILL
  WRITE (*,*) 'GELMAN STORAGE ALLOCATED'
  CALL HRULEFILL


  !$OMP PARALLEL SHARED(DRAWstates,DRAWsvol,GELMANstates,GELMANsvol,DRAWf,GELMANhSigma, DRAWhgap0, GELMANhgap0,DRAWhSigma, DRAWshockslopes,GELMANF)

  !$OMP DO 
  DO k = 1,T
     DO j = 1, Nstates
        call GelmanTest1(GELMANstates(j,k), DRAWstates(j,k,:,:), Nsim, Nstreams)
     END DO
  END DO
  !$OMP END DO

  !$OMP DO 
  DO k = 1,T
     DO j = 1, Nsv
        call GelmanTest1(GELMANsvol(j,k), DRAWsvol(j,k,:,:), Nsim, Nstreams)
     END DO

  END DO
  !$OMP END DO

  call GelmanTest1(GELMANhvarbar(1), DRAWhvarbar(1,:,:), Nsim, Nstreams)


  !$OMP DO 
  DO j = 1, Ngap
     call GelmanTest1(GELMANhgap0(j), DRAWhgap0(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP DO 
  DO j = 1, NhSigma
     call GelmanTest1(GELMANhsigma(j), DRAWhsigma(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP DO 
  DO j = 1, Nf
     call GelmanTest1(GELMANf(j), DRAWf(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP DO 
  DO j = 1, Nshockslopes
     call GelmanTest1(GELMANshockslopes(j), DRAWshockslopes(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP END PARALLEL 

  CALL HRULEFILL
  WRITE (*,*) 'GELMAN STATISTICS ARE DONE!'
  CALL HRULEFILL

  IF (ALL(ABS(GELMANstates - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for the STATES'
  ELSE
     WRITE(*,*) 'STATES: GELMAN FAILURE, Max SRstat=', maxval(GELMANstates)
  END IF
  filename = 'GELMAN.states' // filext
  call savemat(GELMANstates, filename)

  IF (ALL(ABS(GELMANsvol - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for the SVOL'
  ELSE
     WRITE(*,*) 'SVOL: GELMAN FAILURE, Max SRstat=', maxval(GELMANsvol)
  END IF
  filename = 'GELMAN.svol' // filext
  call savemat(GELMANsvol, filename)

  IF (ALL(ABS(GELMANhsigma - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HSIGMA'
  ELSE
     WRITE(*,*) 'HSIGMA: GELMAN FAILURE, Max SRstat=', maxval(GELMANhSigma)
  END IF
  filename = 'GELMAN.hSigma' // filext
  call savevec(GELMANhSigma, filename)


  IF (ALL(ABS(GELMANhgap0 - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HGAP0'
  ELSE
     WRITE(*,*) 'HGAP0: GELMAN FAILURE, Max SRstat=', maxval(GELMANhgap0)
  END IF
  filename = 'GELMAN.hgap0' // filext
  call savevec(GELMANhgap0, filename)

  IF (ALL(ABS(GELMANhvarbar - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HVARBAR'
  ELSE
     WRITE(*,*) 'HVARBAR: GELMAN FAILURE, Max SRstat=', maxval(GELMANhvarbar)
  END IF
  filename = 'GELMAN.hvarbar' // filext
  call savevec(GELMANhvarbar, filename)



  IF (ALL(ABS(GELMANshockslopes - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for SHOCKSLOPES'
  ELSE
     WRITE(*,*) 'SHOCKSLOPES: GELMAN FAILURE, Max SRstat=', maxval(GELMANshockslopes)
  END IF
  filename = 'GELMAN.shockslopes' // filext
  call savevec(GELMANshockslopes, filename)

  IF (ALL(ABS(GELMANf - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for F'
  ELSE
     WRITE(*,*) 'F: GELMAN FAILURE, Max SRstat=', maxval(GELMANf)
  END IF
  filename = 'GELMAN.f' // filext
  call savevec(GELMANf, filename)

  DEALLOCATE (GELMANstates, GELMANsvol, GELMANhgap0, GELMANhSigma, GELMANhvarbar, GELMANshockslopes, GELMANf)

  ! ----------------------------------------------------------------------------
  ! FINISHED: GELMAN
  ! ----------------------------------------------------------------------------

  call hrulefill
  WRITE(*,*) 'DONE, BYE BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  ! ----------------------------------------------------------------------------
  ! CLEANUP
  ! ----------------------------------------------------------------------------

  DEALLOCATE (DRAWmaxlambda, DRAWhSigma, DRAWhgap0, DRAWhvarbar, DRAWstates,DRAWsvol, DRAWshockslopes, DRAWf, STAT=status)
  DEALLOCATE (ydata, yNaN, E0shockslopes, sqrtV0shockslopes, hSigmaT, Eh0, sqrtVh0, minSV, maxSV, SigmaY)

  ! clean out streams
  ! VSLstreams
  errcode = vsldeletestream(VSLdefaultstream)     

  STOP

CONTAINS

  ! -----------------------------------------------------------------
  SUBROUTINE getarguments(datalabel,p,T)

    INTENT(INOUT) datalabel, p, T

    INTEGER :: T, p, counter
    CHARACTER (LEN=100) :: datalabel
    CHARACTER(len=32) :: arg

    counter = 0
    IF (command_argument_count() == 0) THEN
       print *, 'WARNING. No Datalabel specified!'
       print *, 'Using default: ' // datalabel
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, datalabel) 
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') p
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') T
    END IF

  END SUBROUTINE getarguments

  ! -----------------------------------------------------------------
  SUBROUTINE getsettings(datalabel, Ny)

    INTENT(IN) datalabel
    INTENT(OUT) Ny

    INTEGER :: Ny
    CHARACTER (LEN=100) :: datalabel
    CHARACTER (LEN=100) :: filename

    filename = trim(datalabel) // '.settings.txt'
    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

    READ(4,'(t5,i5)') Ny

    CLOSE(UNIT=4)

  END SUBROUTINE getsettings


  ! -----------------------------------------------------------------
  SUBROUTINE readnandata(nanny,filename,Ny,T)
    IMPLICIT NONE

    INTENT(IN) :: filename,T,Ny
    INTENT(INOUT) :: nanny
    CHARACTER (LEN=100) :: filename
    CHARACTER (LEN=500) :: fmtstr

    LOGICAL, DIMENSION(:,:) :: nanny
    INTEGER :: work(Ny)

    INTEGER i, j, T, Ny

    fmtstr = '(I2' // repeat(',I2', Ny-1) // ')'

    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

    DO i=1,T
       READ(4,fmtstr) (work(j), j=1,Ny)
       WHERE (work == 1) 
          nanny(:,i) = .TRUE.
       ELSEWHERE
          nanny(:,i) = .FALSE.
       END WHERE
    END DO

    CLOSE(UNIT=4)

  END SUBROUTINE readnandata
  ! -----------------------------------------------------------------

  ! -----------------------------------------------------------------
  SUBROUTINE readdata(y,filename,Ny,T)
    IMPLICIT NONE

    INTENT(IN) :: filename,T,Ny
    INTENT(INOUT) :: y
    CHARACTER (LEN=100) :: filename
    CHARACTER (LEN=1000) :: fmtstr

    DOUBLE PRECISION, DIMENSION(:,:) :: y
    INTEGER i, T, Ny

    fmtstr = es30d16(Ny)

    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')
    DO i=1,T
       READ(4,fmtstr) y(:,i)
    END DO

    CLOSE(UNIT=4)

  END SUBROUTINE readdata

  ! -----------------------------------------------------------------

END PROGRAM
! -----------------------------------------------------------------

SUBROUTINE thissampler(doDiffuse,T,p,y,yNaN,Ny,DRAWstates,Nstates,Nx,DRAWsvol,Nsv,Ehbar0, Vhbar0,Eh0,sqrtVh0,minSV,maxSV, DRAWhvarbar, hvarbarT, hvarbarDof, DRAWhgap0, DRAWhSigma, NhSigma, hSigmaT, hSigmaDof, DRAWshockslopes, Nshockslopes, E0shockslopes, sqrtV0shockslopes, DRAWf, Nf, Ef0, sqrtVf0, DRAWmaxlambda, Nsim, Burnin, VSLstream,timer)

  use embox, only : hrulefill, savemat, savevec, int2str
  use blaspack, only : eye, vech
  use gibbsbox, only : varmaxroot, igammadraw, drawRW, drawRWcorrelated, iwishcholdraw, bayesVARSV, bayesdiffuseVARSV, stochvolKSC0, SVHcholeskiKSCAR1corplus, SVHdiffusecholeskiKSCAR1, varianceDRAW, vcvcholDrawTR
  use statespacebox, only : DLYAP, samplerA3B3C3nanscalar
  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  INTENT(INOUT) :: DRAWstates,DRAWsvol,DRAWf, DRAWhgap0, DRAWhsigma, DRAWhvarbar, DRAWshockslopes, VSLstream, timer
  INTENT(IN)    :: doDiffuse,T,p,y,yNaN,Ny,Nstates,Nx,Nsv,Ehbar0,Vhbar0,Eh0,sqrtVh0, hvarbarT, hvarbarDof, hSigmaT, hSigmaDof, NhSigma, minSV,maxSV,Nsim,Burnin, Nf, Ef0, sqrtVf0, Nshockslopes, E0shockslopes, sqrtV0shockslopes

  INTEGER :: J, I, K, T, Nsim, Burnin, TotalSim, Nstates, Nx, Ny, Nsv, NhSigma, p, Nbar, Ngap, Nf, status, Nw, Nshockslopes, offsetslopes, these

  logical :: doDiffuse

  type (vsl_stream_state) :: VSLstream
  type(progresstimer) :: timer

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN

  DOUBLE PRECISION, DIMENSION(Nstates,0:T,Nsim) :: DRAWstates
  DOUBLE PRECISION, DIMENSION(Nsv,0:T,Nsim) :: DRAWsvol
  DOUBLE PRECISION, DIMENSION(Nf,Nsim) :: DRAWf
  DOUBLE PRECISION, DIMENSION(NhSigma,Nsim) :: DRAWhSigma
  DOUBLE PRECISION, DIMENSION(Ny,Nsim) :: DRAWhgap0
  DOUBLE PRECISION, DIMENSION(1,Nsim) :: DRAWmaxlambda, DRAWhvarbar
  DOUBLE PRECISION, DIMENSION(Nshockslopes,Nsim) :: DRAWshockslopes

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), A(Nx,Nx,T), Bsv(Nx,1+Ny,T), B(Nx,1+Ny), C(Ny,Nx,T), Ef0(Nf), sqrtVf0(Nf,Nf), iVf0(Nf,Nf), f(Nf), maxlambda, Iy(Ny,Ny), SigmaXstar(Nx,Nx), ygap0variance(Ny * p, Ny * p), gapshock0loadings(Ny * p, Ny)

  DOUBLE PRECISION :: E0shockslopes(Nshockslopes), sqrtV0shockslopes(Nshockslopes,Nshockslopes), iV0shockslopes(Nshockslopes,Nshockslopes), shockslopes(Nshockslopes)
  DOUBLE PRECISION, PARAMETER :: unity = 1.0d0 - 1.0d-6

  DOUBLE PRECISION :: gapdraw(-(p-1):T,Ny), gapshock(1:T,Ny), gapVCV(Ny,Ny,1:T)

  DOUBLE PRECISION :: Ehbar0, Vhbar0, Eh0(Ny), minSV(Nsv), maxSV(Nsv)
  INTEGER :: hSigmaDof, hvarbarDof
  DOUBLE PRECISION, DIMENSION(Ny,Ny) :: sqrtVhshock, sqrtVh0, hSigmaT, hSigma0
  DOUBLE PRECISION :: hvarbarT, hvarbar
  DOUBLE PRECISION, DIMENSION(Nx, 0:T) :: x
  DOUBLE PRECISION, DIMENSION(Nx, T) :: xshock 
  DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: SVol, h
  DOUBLE PRECISION, DIMENSION(Nsv,T)  :: hshock
  DOUBLE PRECISION, DIMENSION(Ny)  :: hgap0
  DOUBLE PRECISION, PARAMETER :: hgaprho = 0.95d0 

  ! Forecasting draws
  DOUBLE PRECISION :: THISh(Nsv), THIShvar(Nsv), THIShgap0(Ny)

  ! stack management
  ! ORG: integer, parameter :: stacksize = 50, maxShakes = 1
  integer, parameter :: stacksize = 10, maxShakes = 1
  double precision, dimension(Nsv,0:T,stacksize) :: PREV_SVOL
  double precision, dimension(Nf,stacksize) :: PREV_F
  double precision, dimension(Ny,Ny,stacksize) :: PREV_HSIGMA
  double precision, dimension(Ny,stacksize) :: PREV_HGAP0
  double precision, dimension(1,stacksize) :: PREV_HVARBAR
  double precision, dimension(Nshockslopes,stacksize) :: PREV_SHOCKSLOPES
  integer, dimension(stacksize) :: stackRegister
  integer :: lastInStack, shakes, stackResetCount
  logical :: OK 
  CHARACTER (LEN=200) :: resetmsg

  ! state space
  INTEGER :: ndxGapCompanionStart, ndxGapCompanionStop

  ! VSL
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  stackRegister = -1 ! dummy value for unregistered stack
  stackResetCount = -1

  Nbar = Ny
  Ngap = Ny
  Nw = 1 + Ny
  TotalSim = Nsim + Burnin
  call eye(Iy)

  ! prepare state space
  A = 0.0d0
  ! unit roots for bar
  FORALL (j=1:Nbar) A(j,j,:) = 1.0d0

  ndxGapCompanionStart = Nbar + 1
  ndxGapCompanionStop  = Nbar + (Ngap * p)

  ! kompanion
  FORALL (j = 1 : Ngap * (p - 1))
     A(Nbar+Ngap+j,Nbar+j,:) = 1.0d0
  END FORALL

  ! C measurement loadings
  C = 0.0d0
  FORALL (j=1:Ny) C(j,j,:) = 1.0d0
  FORALL (j=1:Ny) C(j,Ny+j,:) = 1.0d0
  ! prepare C for missing values
  DO j=1,T
     DO i = 1, Ny
        if (yNaN(i,j)) C(i,:,j) = 0.0d0
        if (yNaN(i,j) .AND. y(i,j) /= 0.0d0 ) then
           write (*,*) 'YNAN PATTERN DOES NOT MATCH ZEROS IN Y'
        end if
     END DO
  END DO

  iVf0 = sqrtVf0
  call DPOTRI('U', Nf, iVf0, Nf, status)
  if (status /= 0) then
     write(*,*) 'DPOTRI ERROR, INFO: ', status, ' [PADDINGTONSAMPLER]'
     stop 1
  end if

  iV0shockslopes = sqrtV0shockslopes
  call DPOTRI('U', Nshockslopes, iV0shockslopes, Nshockslopes, status)
  if (status /= 0) then
     write(*,*) 'DPOTRI ERROR, INFO: ', status, ' [PADDINGTONSAMPLER]'
     stop 1
  end if

  ! prepare prior VCV of states
  Ex0          = 0.0d0
  Ex0(1:Ny)    = 2.0d0
  ! - expressed as lower  triangular-choleski factor (sqrtVx0 * sqrtVx0')
  sqrtVx0 = 0.0d0
  sqrtVx0(1:Nbar,1) = 100.d0                 ! uncertainty about "primary" trend
  FORALL (j=2:Nbar) sqrtVx0(j,j) = 2.0d0     ! uncertainty about other trend given primary





  shakes = 1
  j = 1
  lastInStack = 1 ! redundant, but avoids compiler warning
  resetmsg = 'DRAW 0'

  DO ! while j<= TotalSim

     IF (j == 1) THEN

        stackResetCount = stackResetCount + 1
        call hrulefill
        print *, 'RE-INIT STREAM ', timer%rank, '(reset count = ', stackResetCount, ', ',  trim(resetmsg), ')'
        call initprogressbar(timer,timer%wait,timer%rank)
        call hrulefill


        ! init Trend SV
        DO ! repeat until SV inside bounds
           ! hinno
           if (hvarbarDof == 0) then
              call igammaDraw(hvarbar, .2d0, 2, VSLstream)
           else
              call igammaDraw(hvarbar, hvarbarT, hvarbarDof, VSLstream)
           end if

           ! SV
           CALL drawRW(h(1,:),T,sqrt(hvarbar),Ehbar0,Vhbar0, VSLstream)
           SVol(1,:) = exp(h(1,:) * 0.5d0)
           OK = (minval(SVol(1,:)) .ge. minSV(1)) .AND. (maxval(SVol(1,:)) .le. maxSV(1))
           IF (OK) EXIT
        END DO

        ! - hSigma
        if (hSigmaDof > 0) then
           call iwishcholDraw(sqrtVhshock, hSigmaT, hSigmaDof, Nsv-1, VSLstream)
        else
           call eye(hSigma0, 0.1d0 ** 2)
           call iwishcholDraw(sqrtVhshock, hSigma0, Nsv - 1 + 2, Nsv-1, VSLstream)
        end if

        ! draw SV
        CALL drawRWcorrelated (h(2:Nsv,:),Nsv-1,T,sqrtVhshock,Eh0,sqrtVh0,VSLstream)
        SVol(2:Nsv,:) = exp(h(2:Nsv,:) * 0.5d0)
        WHERE (SVol .gt. 10.0d0) SVol = 10.0d0
        WHERE (SVol .lt. 1.0d-2) SVol = 1.0d-2



        ! init VAR coefficients f
        ! print *, 'STREAM ', timer%rank, 'starting init VAR f ... '
        DO
           errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nf, f, 0.0d0, 1.0d0)
           CALL DTRMV('U', 'T', 'N', Nf, sqrtVf0, Nf, f, 1)
           f = f + Ef0
           call VARmaxroot(maxlambda, f, Ngap, p)
           OK = maxlambda < unity
           IF (OK) EXIT

        END DO
        ! print *, 'STREAM ', timer%rank, 'done init VAR f'


        ! init shockslopes
        errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nshockslopes, shockslopes, 0.0d0, 1.0d0)
        CALL DTRMV('U', 'T', 'N', Nshockslopes, sqrtV0shockslopes, Nshockslopes, shockslopes, 1)
        shockslopes = shockslopes + E0shockslopes


        ! init stack
        lastInStack = 1

        PREV_HVARBAR(:,lastInStack)     = hvarbar
        PREV_SVOL(:,:,lastInStack)      = SVol
        PREV_HSIGMA(:,:,lastInStack)    = sqrtVhshock
        PREV_HGAP0(:,lastInStack)       = Eh0
        PREV_F(:,lastInStack)           = f
        PREV_SHOCKSLOPES(:,lastInStack) = shockslopes

        stackRegister = -1
        stackRegister(lastInStack)      = j

     ELSE ! j > 1

        hvarbar     = PREV_HVARBAR(1,lastInStack)   
        SVol        = PREV_SVOL(:,:,lastInStack)
        sqrtVhshock = PREV_HSIGMA(:,:,lastInStack)      
        hgap0       = PREV_HGAP0(:,lastInStack)      
        f           = PREV_F(:,lastInStack) 
        shockslopes = PREV_SHOCKSLOPES(:,lastInStack)


     END IF ! j == 1

     ! check inputs
     if (any(SVol <= 0.0d0)) then
        WRITE(*,*) 'negative SVOL?!'
        STOP 1
     end if


     ! ---------------------------------------------------------------------------
     ! SMOOTHING SAMPLER
     ! ---------------------------------------------------------------------------

     ! a) VAR coefficients
     ! Fill in VAR coefficients
     FORALL (i=1:T) A(Nbar+1:Nbar+Ngap,Nbar+1:Nbar+Ngap*p,i) = transpose(reshape(f, (/ Ngap * p, Ngap /)))


     ! b1) shockslopes 
     ! note: Since B gets multiplied by SVol at each iteration, need to build B from scratch for every iteration j
     B = 0.0d0
     B(1:Nbar,1) = 1.0d0
     call eye(B(Nbar+1:Nbar+Ngap,2:1+Ngap))
     offsetslopes = 0
     DO i = 2, Ngap 
        these = i-1 
        ! slopes in row i have index offsetslopes+1:offsetslopes + these
        B(Nbar+i,2:1+these) = shockslopes(offsetslopes+1:offsetslopes+these)
        offsetslopes = offsetslopes + these
     END DO

     ! c) SV (scaling columns of B)
     ! note: could also call SV sampler
     Bsv = 0.0d0
     FORALL (k=1:T,i=1:Nsv) Bsv(:,i,k) = B(:,i) * SVol(i,k)

     ! update prior Variance for initial states 
     FORALL (i=1:Ngap) gapshock0loadings(:,i) = B(ndxGapCompanionStart:ndxGapCompanionStop,1+i) * SVol(1+i,0)
     CALL DLYAP(ygap0variance, A(ndxGapCompanionStart:ndxGapCompanionStop,ndxGapCompanionStart:ndxGapCompanionStop,1), gapshock0loadings, Ngap * p, Ngap, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DLYAP error (ygap0variance)', errcode
        stop 1
     end if

     ! Factorize the unconditional variance
     CALL DPOTRF('L', Ngap * p, ygap0variance, Ngap * p, errcode)
     if (errcode /= 0) then
        write (*,*) 'DPOTRF error (ygap0variance)', errcode
        stop 1
     end if
     ! zero out the upper triangular
     FORALL (i=2:Ngap*p) ygap0variance(1:i-1,i) = 0.0d0
     ! fill it in
     sqrtVx0(ndxGapCompanionStart:ndxGapCompanionStop,ndxGapCompanionStart:ndxGapCompanionStop) = ygap0variance


     ! call smoothing sampler
     CALL samplerA3B3C3NaNscalar(x,xshock,SigmaXstar,y,yNaN,T,Ny,Nx,Nw,A,Bsv,C,Ex0,sqrtVx0,VSLstream,status)
     ! SigmaXstar is a dummy here 

     if (status /= 0) then
        print *, 'status is', status
        OK = .FALSE.

        ! write out debug data
        call savemat(A(:,:,1), 'A.debug')
        call savemat(Bsv(:,:,1), 'Bsv.debug')
        call savemat(B, 'B.debug')
        call savemat(SVol, 'Svol.debug')

        print *, j
        stop 33

     else
        OK = .TRUE.

        ! fill in forecasting draws
        THISh                         = 2* log(SVol(:,T))
        THIShvar(1)                   = hvarbar
        forall (i=1:Ny) THIShvar(1+i) = sqrtVhshock(i,i) ** 2
        THIShgap0                     = hgap0


     end if


     ! if (any(isnan(x))) then
     !    write(*,*) 'NaN draws from Kalman Sampler'
     !    stop 1
     ! end if

     IF (OK) THEN
        ! ---------------------------------------------------------------------------
        ! VAR Step
        ! ---------------------------------------------------------------------------

        ! a1: construct gapVCV
        gapVCV = 0.0d0
        DO i = 1,T
           call dsyrk('U', 'n', Ngap, Nw, 1.0d0, Bsv(Nbar+1:Nbar+Ngap,:,i), Ngap, 0.0d0, gapVCV(:,:,i), Ngap)
        END DO

        ! a2: invert gapVCV
        DO i = 1,T
           call DPOTRF('U', Ngap, gapVCV(:,:,i), Ngap, status)
           if (status /= 0) then
              call savemat(Bsv(Nbar+1:Nbar+Ngap,:,i), 'Bsv.debug')
              call savemat(B, 'B.debug')

              ! OK = ALL(SVol(4,:) >= minSV(4)) .AND. ALL(SVol(4,:) <= maxSV(4))
              ! print *, 'OK', OK

              print *, 'MCMC step', j
              print *, 't', i
              print *, 'Svol(t)', SVol(:,i)
              write(*,*) 'DPOTRF ERROR, INFO: ', status, ' [GAPVCV - PADDINGTONSAMPLER]'
              stop 1

           end if
           call DPOTRI('U', Ngap, gapVCV(:,:,i), Ngap, status)
           if (status /= 0) then
              write(*,*) 'DPOTRI ERROR, INFO: ', status, ' [GAPVCV - PADDINGTONSAMPLER]'
              stop 1
           end if
        END DO


        ! b) collect lagged gaps (note: will be zero under tight prior)
        FORALL (i=1:p) gapdraw(1-i,:) = x(Nbar+(i-1)*Ngap+1:Nbar+i*Ngap,0)
        gapdraw(1:T,:)                = transpose(x(Nbar+1:Nbar+Ngap,1:T))


        ! d) draw VAR coefficients
        ! recall: gapVCV contains inverse of gapVCV at this point
        if (doDiffuse) then
           call bayesdiffuseVARSV(f, gapshock, p, gapdraw, Ngap, T, gapVCV, VSLstream)
        else
           call bayesVARSV(f, gapshock, p, gapdraw, Ngap, T, gapVCV, Ef0, iVf0, VSLstream)
        end if


        if (any(isnan(f))) then
           write(*,*) 'NaN draws from bayesVAR'
           stop 1
        end if

        ! d) check stability
        call VARmaxroot(maxlambda, f, Ngap, p)

        ! print *, 'maxroot', maxlambda
        ! stop 1

        OK = maxlambda < unity
        IF (.NOT. OK) resetmsg = 'VAR'

        xshock(Nbar+1:Nbar+Ngap,:) = transpose(gapshock)

        ! ---------------------------------------------------------------------------
     END IF

     h = 2.0d0 * log(SVol)

     IF (OK) THEN

        ! ---------------------------------------------------------------------------
        ! Trend SV
        ! ---------------------------------------------------------------------------
        CALL stochvolKSC0(h(1,:),  xshock(1,:), T, sqrt(hvarbar), Ehbar0, Vhbar0, VSLstream)
        hshock(1,:) = h(1,1:T) - h(1,0:T-1)

        SVol(1,:) = exp(h(1,:) * 0.5d0)

        OK = ALL(SVol(1,:) >= minSV(1)) .AND. ALL(SVol(1,:) <= maxSV(1))
        IF (.NOT. OK) resetmsg = 'SVOL 1'




     END IF


     IF (OK) THEN

        ! ---------------------------------------------------------------------------
        ! Gap SV block
        ! ---------------------------------------------------------------------------

        if (doDiffuse) then
           call SVHdiffusecholeskiKSCAR1(T, Ny, SVol(2:Nsv,:), h(2:Nsv,:), hshock(2:Nsv,:), hgap0, hgaprho, shockslopes, Nshockslopes, xshock(Nbar+1:Nbar+Ngap,:), sqrtVhshock, Eh0, sqrtVh0, VSLstream)
        else
           call SVHcholeskiKSCAR1corplus(T, Ny, SVol(2:Nsv,:), h(2:Nsv,:), hshock(2:Nsv,:), hgap0, hgaprho, shockslopes, Nshockslopes, xshock(Nbar+1:Nbar+Ngap,:), E0shockslopes, iV0shockslopes, sqrtVhshock, Eh0, sqrtVh0, VSLstream)
        end if

        OK = .true.

        DO i=2,Nsv
           IF (ANY(SVol(i,:) < minSV(i))) THEN
              OK = .FALSE.
              resetmsg = 'SVOL Gap ' // trim(int2str(i)) // ' too low'
              EXIT
           END IF
           IF (ANY(SVol(i,:) > maxSV(i))) THEN
              OK = .FALSE.
              resetmsg = 'SVOL Gap ' // trim(int2str(i)) // ' too high'
              EXIT
           END IF
        END DO

     END IF



     IF (OK) THEN ! move forward


        ! ---------------------------------------------------------------------------
        ! vol of Gap h
        ! ---------------------------------------------------------------------------
        call varianceDraw(hvarbar, hvarbarT, hvarbarDof, hshock(1,:), T, VSLstream)
        call vcvcholDrawTR(sqrtVhshock, hSigmaT, hSigmaDof, transpose(hshock(2:Nsv,:)), T, Nsv-1, VSLstream)

        ! a) update stack
        lastInStack = minloc(stackRegister,1) ! find a free element on the stack
        !        WRITE(*,*) "setting lastInStack =", lastInStack
        stackRegister(lastInStack) = j + 1 ! notice the +1
        PREV_SVOL(:,:,lastInStack)      = SVol
        PREV_HSIGMA(:,:,lastInStack)    = sqrtVhshock
        PREV_HVARBAR(:,lastInStack)     = hvarbar
        PREV_HGAP0(:,lastInStack)       = hgap0
        PREV_F(:,lastInStack)           = f
        PREV_SHOCKSLOPES(:,lastInStack) = shockslopes

        ! b) store Draws after Burnin
        IF (j > BURNIN) THEN
           DRAWstates(1:Ny,:,j-BURNIN) = x(1:Nbar,:)
           DRAWstates(Ny+1:Nstates,:,j-BURNIN) = x(Nbar+1:Nbar+Ngap,:)

           DRAWf(:,j-BURNIN)           = f
           DRAWshockslopes(:,j-BURNIN) = shockslopes
           DRAWmaxlambda(:,j-BURNIN)   = maxlambda
           DRAWsvol(:,:,j-BURNIN)      = SVol
           DRAWhvarbar(:,j-BURNIN)      = hvarbar
           DRAWhgap0(:,j-BURNIN)      = hgap0
           call vech(DRAWhsigma(:,j-BURNIN), sqrtVhshock)

           ! CONSTRUCT FORECASTS AND STORE THEM

           ! fill in variables as per the Kalman step

        END IF
        j = j + 1

     ELSE ! NOT(OK): shake or move back

        OK = .TRUE.

        IF (shakes < maxshakes) THEN
           shakes = shakes + 1
           WRITE (*,'("Stream ", i2, ": shaking at step ", i5)') timer%rank, j

        ELSE ! move back

           shakes = 1
           IF (j > 1) THEN
              ! WRITE (*,'("Stream ", i2, ": moving back at step ", i5, " -- " a20)') timer%rank, j, resetmsg
              stackRegister(lastInStack) = -1
              lastInStack = maxloc(stackRegister,1)

              ! WRITE(*,*) "stackRegister=", stackRegister(lastInStack), 'j=', j
              j = j - 1
              if (stackRegister(lastInStack) < 1) then
                 j = 1
              else
                 if (stackRegister(lastInStack) .ne. j) then
                    print *, 'stack problem', 'j=', j, 'stackregister=', stackregister(lastinstack)
                 end if
              end if



              ! IF (lastInStack < 1) THEN 
              !    WRITE(*,*) "lastInStack below 1 -- aborting ..."
              !    CALL EXIT(1)
              ! END IF

              ! ELSE ! j = 1
              !    stackregister = -1
              !    IF (any(stackregister .ge. 0.0d0)) THEN
              !       call savevec(dble(stackregister), 'stackregister.debug')
              !       print *, 'stack problem'
              !       stop 1
              !    end if

           END IF ! j > 1

        END IF ! shakes 

     END IF ! OK

     if (j > TotalSim) EXIT

     CALL progressbar(dble(j) / dble(TotalSim), timer)


  END DO ! while j < TotalSim


  print *, 'Stream ', timer%rank, ' had ', stackResetCount, 'resets'

END SUBROUTINE thissampler
