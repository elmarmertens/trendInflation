PROGRAM main

  ! particle filter for commontrend model (no trend slopes)
  ! packed storage for SigmaXposterior in particlefilter

  USE embox, only : hrulefill, loadmatrix, loadarray1, savemat, savematlogical, savevec, storeEstimates, loft, timestampstr, es30d16, int2str
  USE blaspack, only : vech, ivech
  USE gibbsbox, only : drawNDXpdf, drawNDXsysresample

  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  INTEGER, PARAMETER :: parameterColumns = 12, meanCol = 1 ! meanCol is used to pick parameters from input file

  logical, parameter :: doTimestamp = .false., doSmoother = .false., doGains = .true.

  INTEGER :: p, Ny, Nx, Nw, Nsv, Nf, NhgapSigma, Nstates, Nshockslopes

  INTEGER :: Nparticles, Nsmoother, NsmootherX, Nmixturedraws, Ndraws
  INTEGER :: T,i,j,k,status ! ,q

  ! filter particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PARTICLEweights, DRAWllf, DRAWlike
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWsvol, DRAWxhat, DRAWxsig
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWxgain
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XBAR

  ! smoother particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: SMOOTHERsvol
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: SMOOTHERx

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sqrtVx0, sqrtVhgap
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: f, shockslopes, Ex0, hgap0
  DOUBLE PRECISION, PARAMETER :: hgaprho = 0.95d0 
  DOUBLE PRECISION :: hbarinno, Ehbar0, Vhbar0

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: loglike, loglike2
  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Xdraws
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: theta2
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: theta1

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ndx

  TYPE(progresstimer) :: timer ! , stopwatch
  CHARACTER (LEN=200) :: filename, datafile, nandatafile, fileXT, datalabel, parameterlabel, parameterXT

  ! VSL Random Stuff"
  type (vsl_stream_state) :: VSLstream
  integer :: seed
  integer :: brng
  integer :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OPEN MP
  INTEGER :: NTHREADS !, TID

  ! ----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  ! ----------------------------------------------------------------------------------
  ! runtime parameters :start:
  ! first: set default values

  ! thorough
  Nparticles    = 4 * (10 ** 4)
  ! Nparticles    = 10 ** 5
  Nmixturedraws = 10 ** 4
  Nsmoother     = 10 ** 4
  NsmootherX    = 100

 
  ! quick
  Nparticles    = 10 ** 3
  Nmixturedraws = 10 ** 3
  Nsmoother     = 10 ** 3
  NsmootherX    = 10
  
  datalabel        = 'INFTRM'
  parameterlabel   = 'DEFAULT'


  T = 0
  call getarguments(datalabel, T, parameterlabel) 
  IF (parameterlabel == 'DEFAULT')  parameterlabel = datalabel
  call getsettings(parameterlabel,Ny)


  p      = 12
  
  Nshockslopes = Ny * (Ny - 1) / 2
  NhgapSigma   = Ny * (Ny + 1) / 2
  Nstates      = Ny 
  Nx           = Ny + Ny * p
  Nsv          = 1 + Ny
  Nw           = 1 + Ny
  Nf           = p * Ny * Ny


  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  call hrulefill
  if (doSmoother) then
     print *, 'Particle Filter estimation of PaddingtonGAPSV'
  else
     print *, 'Particle Filter (and Smoother) estimation of Paddington'
  end if
  

  ! INIT OMP
  NTHREADS = 1
  !$OMP PARALLEL SHARED(NTHREADS)
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  print *, "Number of Threads:", NTHREADS

  ! VSL
  brng    = vsl_brng_mt19937
  seed    = 0
  errcode = vslnewstream(VSLstream, vsl_brng_mt2203, seed)  
  
  ! runtime parameters :end: 

  ! CONSTRUCT FILE EXTENTSIONS
  fileXT = '.particles.' // trim(datalabel) // '.gapSV.dat'
  if (doTimeStamp) filext = '.' // timestampstr() //  filext

  datafile    = trim(datalabel) // '.yData.txt'
  nandatafile = trim(datalabel) // '.yNaN.txt'

  ! read data
  if (T == 0) then
     T = loft(datafile) 
     IF (T < 10) THEN
        print *, 'Less than 10 observations in input file!', datafile
        STOP 1
     END IF
  end if


  parameterXT   = '.notrendslopes.' // trim(parameterlabel) // '.T' // trim(int2str(T))  // '.gapSV.dat'

  ALLOCATE (y(Ny,T), yNaN(Ny,T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Y)'
  END IF

  ! print *, 'trying to read', T, 'obs from', datafile
  CALL readdata(y,datafile,Ny,T)
  CALL readnandata(yNaN,nandatafile,Ny,T)

  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data= ' // datalabel
  print *, 'parameters= ' // parameterlabel
  print *, 'Ny= ', Ny
  print *, 'T= ', T
  print *, 'Nparticles= ', Nparticles
  print *, 'Nsmoother= ', Nsmoother
  print *, 'p= ', p
  CALL HRULEFILL

  ! Model parameters and priors
  ALLOCATE (hgap0(Ny), sqrtVhgap(Ny,Ny))
  ALLOCATE (f(Nf), shockslopes(Nshockslopes))

  ! SV bar moments
  Vhbar0   = 1.0d0
  Ehbar0   = -2.0d0 * log(12.0d0) - Vhbar0 * 0.25d0


  ! Linear prior
  ALLOCATE (Ex0(Nx), sqrtVx0(Nx,Nx))
  Ex0       = 0.0d0
  Ex0(1:Ny) = 2.0d0

 ! sqrtVx0, expressed as lower  triangular-choleski factor (sqrtVx0 * sqrtVx0')
 sqrtVx0 = 0.0d0
 sqrtVx0(1:Ny,1) = 100.d0                 ! uncertainty about "primary" trend
 FORALL (j=2:Ny) sqrtVx0(j,j) = 2.0d0     ! uncertainty about other trend given primary

 


  ! ---------------------------------
  ! LOAD PARAMETERS
  ! ---------------------------------

  ! 1) parameter files with quantiles

  ! f
  allocate(theta2(Nf, parameterColumns))
  filename = 'F' // trim(adjustl(parameterXT))
  call loadmatrix(theta2, filename, Nf, parameterColumns)
  f = theta2(:,meanCol)
  deallocate(theta2)

  print *, 'f'
  print *, f
  call hrulefill


  ! shockslopes
  filename = 'SHOCKSLOPES' // trim(adjustl(parameterXT))
  allocate(theta2(Nshockslopes, parameterColumns))
  call loadmatrix(theta2, filename, Nshockslopes, parameterColumns)
  shockslopes = theta2(:,meanCol) 
  deallocate(theta2)
  print *, 'shockslopes'
  print *, shockslopes
  call hrulefill

  ! sqrtVhgap
  sqrtVhgap  = 0.0d0
  filename = 'HSIGMA' // trim(adjustl(parameterXT))
  allocate(theta2(NhgapSigma,parameterColumns))
  call loadmatrix(theta2, filename, NhgapSigma, parameterColumns)
  call ivech(sqrtVhgap, theta2(:,meanCol))
  ! recall: sqrtVhgap is upper-triangular draw from inverse wishart, i.e. Vhgap = sqrtVhgap * sqrtVhgap'
  ! call savemat(sqrtVhgap, 'sqrtVh0.debug')
  deallocate(theta2)
  print *, 'sqrtVhgap'
  print *, sqrtVhgap
  call hrulefill




  ! 2) parameters stored as draws



  ! hbarinno
  filename = 'HVARBAR' // trim(adjustl(parameterXT))
  Ndraws = loft(filename)
  allocate(theta1(Ndraws))
  call loadarray1(theta1, filename, Ndraws)
  if (meanCol == 1) then
     hbarinno = sum(sqrt(theta1)) / dble(Ndraws)
  else ! do median
     CALL dlasrt('I', Ndraws, theta1, status)
     hbarinno = sqrt(theta1(floor(real(Ndraws) * 0.5)))
  end if
  deallocate(theta1)
  print *, 'hbarinno'
  print *, hbarinno
  call hrulefill


  ! hgap0
  filename = 'HGAP0' // trim(adjustl(parameterXT))
  Ndraws = loft(filename)
  allocate(theta2(Ndraws, Ny))
  call loadmatrix(theta2, filename, Ndraws, Ny)
  if (meanCol == 1) then
     hgap0 = sum(theta2, 1) / dble(Ndraws)
  else ! do median
     do j=1,Ny
     CALL dlasrt('I', Ndraws, theta2(:,j), status)
     end do
     hgap0 = theta2(floor(real(Ndraws) * 0.5),:)
  end if
  deallocate(theta2)
  print *, 'hgap0'
  print *, hgap0
  call hrulefill

  ! allocate memory for draws
  ALLOCATE (PARTICLEweights(0:T,Nparticles),DRAWllf(T,Nparticles), DRAWlike(T,Nparticles), DRAWxhat(Nstates,T,Nparticles), DRAWxsig(Nstates,T,Nparticles), DRAWsvol(Nsv,0:T,Nparticles), DRAWxgain(Nstates,Ny,T,Nparticles), XBAR(Nstates,T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (draws)'
  END IF

  PARTICLEweights = 0.0d0
  DRAWllf         = 0.0d0
  DRAWlike        = 0.0d0
  DRAWxhat        = 0.0d0
  DRAWxsig        = 0.0d0
  DRAWsvol        = 0.0d0
  DRAWxgain       = 0.0d0
  XBAR            = 0.0d0

  CALL particlefilter(T, Ny, y, yNaN, Nparticles, XBAR, PARTICLEweights, DRAWllf, DRAWlike,  DRAWxhat, DRAWxsig, DRAWxgain, Nstates, Nx, Nw, Ex0, sqrtVx0, DRAWsvol, Nsv, Ehbar0, Vhbar0, hbarInno, hgap0, hgaprho, sqrtVhgap, f, p, shockslopes, Nshockslopes, VSLstream,timer)


  CALL HRULEFILL
  WRITE (*,*) 'PARTICLE FILTER IS DONE!'
  CALL HRULEFILL

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a20)') 'TIME: ', timestampstr()
  WRITE(4,'(a20,a20)') 'Data: ', datalabel
  WRITE(4,'(a20,a20)') 'Parameters: ', parameterlabel
  WRITE(4,'(a40)') repeat('-',40)
  WRITE(4,'(a20,I20)') 'Nparticles: ', Nparticles
  WRITE(4,'(a20,I20)') 'p: ', p
  CLOSE(UNIT=4)
  CALL HRULEFILL

  ! ----------------------------------------------------------------------------
  ! STORE
  ! ----------------------------------------------------------------------------

  CALL HRULEFILL
  WRITE (*,*) 'STARTING W/STORAGE !!!'
  CALL HRULEFILL

  ! debugging output
  ! call savemat(DRAWllf, 'DRAWllf.debug')
  ! call savemat(DRAWsvol(1,:,:), 'DRAWsvol.debug')
  ! call savemat(DRAWxhat(:,:,1), 'DRAWxhat.debug')
  ! call savemat(PARTICLEweights, 'PARTICLEweights.debug')

  ! STORE ESTIMATES
  ! Note: manual reshape avoids segmentation faults

  filename = 'YDATA' // filext
  call savemat(y, filename)

  filename = 'YNAN' // filext
  call savematlogical(yNaN, filename)


  ! LIKELIHOOD
  filename = 'LOGLIKE' // filext
  ALLOCATE (loglike(T), loglike2(T), STAT=status)
  loglike  = log(sum(exp(DRAWllf),2) / Nparticles)
  call savevec(loglike, filename)
  call hrulefill
  WRITE (*,*) 'STORED LFF'
  WRITE (*,*) '... the loglikelihood is ', sum(loglike)
  WRITE (*,*) '... w/average contribution ', sum(loglike) / T

  filename = 'LOGLIKE2' // filext
  loglike2  = log(sum(exp(DRAWlike),2) / Nparticles)
  call savevec(loglike2, filename)
  WRITE (*,*) 'STORED LIKE'
  WRITE (*,*) '... the LIKE is ', sum(loglike2)
  WRITE (*,*) '... w/average contribution to the LIKE ', sum(loglike2) / T
  call hrulefill
  DEALLOCATE (DRAWllf,DRAWlike,loglike2)


  ALLOCATE (theta1(T), STAT=status)

  ! ESS
  filename  = 'ESS' // filext
  theta1    = 1 / sum(PARTICLEweights(1:T,:) ** 2, 2) / Nparticles 
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ESS'


  ! store trends
  DO i=1,Ny
     filename  = 'TAUHAT' // trim(int2str(i)) // filext
     theta1 = sum(PARTICLEweights(1:T,:) * DRAWxhat(i,:,:), 2)
     call savevec(theta1, filename)
     WRITE (*,*) 'STORED TAUHAT', i
  END DO

  ! ! store gaps
  ! DO i=1,Ny
  !    filename  = 'GAPHAT' // trim(int2str(i)) // filext
  !    theta1 = sum(PARTICLEweights(1:T,:) * DRAWxhat(Ny+i,:,:), 2)
  !    call savevec(theta1, filename)
  !    WRITE (*,*) 'STORED GAPHAT', i
  ! END DO

  DEALLOCATE (theta1)

  ! draw distribution for linear states
  ALLOCATE (ndx(Nmixturedraws,T),Xdraws(Nstates,T,Nmixturedraws))

  ! sample particle indices (note: these will also be used later for the SV particles)
  print *, 'Drawing particles indices  ...'
  DO j=1,T
     call drawNDXpdf(ndx(:,j), Nmixturedraws, PARTICLEweights(j,:), Nparticles, VSLstream)
  END DO
  print *, 'Done drawing particles indices.'

  print *, 'Drawing Xdraws normals ...'
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Nstates * Nmixturedraws, Xdraws, 0.0d0, 1.0d0)
  print *, 'done drawing Xdraws normals.'
  
! j = 1
! call savevec(dble(ndx(:,j)), 'ndx.debug')
! call savemat(DRAWxhat(:,j,:), 'DRAWxhat.debug')
! call savemat(DRAWxsig(:,j,:), 'DRAWxsig.debug')
! call savemat(Xdraws(:,j,:), 'Zdraws.debug')

  FORALL (i=1:Nstates,j=1:T,k=1:Nmixturedraws) Xdraws(i,j,k) = DRAWxhat(i,j,ndx(k,j)) + DRAWxsig(i,j,ndx(k,j)) * Xdraws(i,j,k) 

! j = 1
! call savemat(Xdraws(:,j,:), 'Xdraws.debug')
! stop 11

  ! store trend draws
  DO i=1,Ny
     filename  = 'TAU' // trim(int2str(i)) // filext
     CALL storeEstimates(Xdraws(i,:,:),T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED TAU', i
  END DO

  ! ! store gap draws
  ! DO i=1,Ny
  !    filename  = 'GAP' // trim(int2str(i)) // filext
  !    CALL storeEstimates(Xdraws(Ny+i,:,:),T,Nmixturedraws,filename)
  !    WRITE (*,*) 'STORED GAP', i
  ! END DO

  
  filename  = 'XBAR' // filext
  CALL savemat(XBAR, filename)
  WRITE (*,*) 'STORED XBAR'
  

  DEALLOCATE (Xdraws, DRAWxhat, DRAWxsig, XBAR)



  ! 2D Gain Matrices
  if (doGains) then
     ! ALLOCATE (theta2(T,Nmixturedraws))

     ! ! store trend gains 
     ! DO i=1,Ny
     !    DO q=1,Ny
     !       filename  = 'GAIN' // trim(int2str(q)) // 'TAU' // trim(int2str(i)) // filext
     !       FORALL (j=1:T,k=1:Nmixturedraws) theta2(j,k) = DRAWxgain(i,q,j,ndx(k,j)) 
     !       CALL storeEstimates(theta2,T,Nmixturedraws,filename)
     !       WRITE (*,*) 'STORED GAIN OF Y', q, 'ON TAU', i
     !    END DO
     ! END DO
     ! ! store gap gains 
     ! DO i=1,Ny
     !    DO q=1,Ny
     !       filename  = 'GAIN' // trim(int2str(q)) // 'GAP' // trim(int2str(i)) // filext
     !       FORALL (j=1:T,k=1:Nmixturedraws) theta2(j,k) = DRAWxgain(Ny+i,q,j,ndx(k,j)) 
     !       CALL storeEstimates(theta2,T,Nmixturedraws,filename)
     !       WRITE (*,*) 'STORED GAIN OF Y', q, 'ON GAP', i
     !    END DO
     ! END DO


     ! DEALLOCATE (theta2)

     ! store analytical moments of gain
     ALLOCATE (theta2(T,Ny))
     ! trend gains 
     DO i=1,Ny
        filename  = 'GAINTAUHAT' // trim(int2str(i)) // filext
        FORALL (j=1:Ny) theta2(:,j) = sum(PARTICLEweights(1:T,:) * DRAWxgain(i,j,:,:), 2)
        call savemat(theta2, filename)
        WRITE (*,*) 'STORED GAINTAUHAT', i
     END DO
     ! ! gap gains 
     ! DO i=1,Ny
     !    filename  = 'GAINGAPHAT' // trim(int2str(i)) // filext
     !    FORALL (j=1:Ny) theta2(:,j) = sum(PARTICLEweights(1:T,:) * DRAWxgain(Ny+i,j,:,:), 2)
     !    call savemat(theta2, filename)
     !    WRITE (*,*) 'STORED GAINGAPHAT', i
     ! END DO
     DEALLOCATE (theta2)
  end if ! doGains
  DEALLOCATE(DRAWxgain)


  ! 2) Nparticle draws for the other particles
  ALLOCATE (theta1(T), STAT=status)
  DO i=1,Nsv
     filename  = 'SVHAT' // trim(int2str(i)) // filext
     theta1 = sum(PARTICLEweights(1:T,:) * DRAWsvol(i,1:T,:), 2)
     call savevec(theta1, filename)
     ! CALL storeEstimates(DRAWsvol(i,1:T,:),T,Nparticles,filename)
     WRITE (*,*) 'STORED SVHAT', i
  END DO
  DEALLOCATE(theta1)

  ! draw distribution
  ALLOCATE (theta2(T,Nmixturedraws))
  DO i=1,Nsv
     filename  = 'SV' // trim(int2str(i)) // filext
     FORALL (j=1:T,k=1:Nmixturedraws) theta2(j,k) = DRAWsvol(i,j,ndx(k,j)) 
     CALL storeEstimates(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED SV', i
  END DO

  DEALLOCATE(theta2)
  DEALLOCATE(ndx)

  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE FILTER
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  ! CLEANUP FILTER
  ! ----------------------------------------------------------------------------

  call hrulefill
  WRITE (*,*) 'LFF:'
  WRITE (*,*) '... the loglikelihood is ', sum(loglike)
  WRITE (*,*) '... w/average contribution ', sum(loglike) / T
  call hrulefill
  DEALLOCATE (loglike)

  IF (.NOT. doSmoother) THEN

     DEALLOCATE (PARTICLEweights, DRAWsvol)
     DEALLOCATE (y, yNaN)
     DEALLOCATE (hgap0, sqrtVhgap)
     DEALLOCATE (f, shockslopes)
     DEALLOCATE (Ex0, sqrtVx0)

     ! VSLstreams
     errcode = vsldeletestream(VSLstream)     
     call hrulefill
     WRITE(*,*) 'DOING ONLY FILTERING.'
     WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
     call hrulefill

     stop 0

  END IF


  ! ----------------------------------------------------------------------------
  ! SMOOTHER
  ! ----------------------------------------------------------------------------
  ALLOCATE (SMOOTHERsvol(Nsv,0:T,Nsmoother))
  
  call hrulefill
  print *, 'STARTING w/SMOOTHER ...'
  CALL particlesmoother(T, Nsmoother, Nsv, SMOOTHERsvol, Nparticles, PARTICLEweights, DRAWsvol, hbarinno, Ny, hgap0, hgaprho, sqrtVhgap, VSLstream,timer)
  print *, 'STARTING w/SMOOTHER STORAGE ...'
  DO i=1,Nsv
     filename  = 'smootherSV' // trim(int2str(i)) // filext
     CALL storeEstimates(SMOOTHERsvol(i,1:T,:),T,Nsmoother,filename)
     WRITE (*,*) 'SMOOTHER: STORED SV', i
  END DO
  DEALLOCATE (PARTICLEweights, DRAWsvol)
  DEALLOCATE (hgap0, sqrtVhgap)

  print *, '... Smoother X ...'
  ALLOCATE (SMOOTHERx(Nstates,0:T,NsmootherX,Nsmoother))
  SMOOTHERx = 0.0d0
  call particleSmootherX(T, Ny, y, yNaN, Nsv, Nsmoother, SMOOTHERsvol, NsmootherX, SMOOTHERx, Nstates, Nx, Nw, Ex0, sqrtVx0, f, p, shockslopes, Nshockslopes)
  DEALLOCATE (f, shockslopes)
  DEALLOCATE (SMOOTHERsvol)


  ! STORE SMOOTHER X
  ALLOCATE (theta2(T,Nsmoother * NsmootherX))
  theta2 = 0.0d0


  ! store trends
  DO i=1,Ny
     filename  = 'smootherTAU' // trim(int2str(i)) // filext
     FORALL (j=1:NsmootherX,k=1:Nsmoother) theta2(:,(k-1) * NsmootherX + j) = SMOOTHERx(i,1:T,j,k)
     CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
     WRITE (*,*) 'SMOOTHER STORED TAU', i
  END DO

  ! ! store gaps
  ! DO i=1,Ny
  !    filename  = 'smootherGAP' // trim(int2str(i)) // filext
  !    FORALL (j=1:NsmootherX,k=1:Nsmoother) theta2(:,(k-1) * NsmootherX + j) = SMOOTHERx(Ny+i,1:T,j,k)
  !    CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
  !    WRITE (*,*) 'SMOOTHER: STORED GAP', i
  ! END DO

  ! FINISHED STORING SMOOTHER X
  DEALLOCATE (SMOOTHERx, theta2)

  print *, '... DONE w/SMOOTHER.'
  call hrulefill

  ! ----------------------------------------------------------------------------
  ! FINAL CLEANUP
  ! ----------------------------------------------------------------------------

  DEALLOCATE (y, yNaN)
  DEALLOCATE (Ex0, sqrtVx0)

  ! VSLstreams
  errcode = vsldeletestream(VSLstream)     
  ! DO j = 0, NTHREADS-1
  !    errcode = vsldeletestream(VSLstreams(j))     
  ! END DO
  ! DEALLOCATE (VSLstreams)

  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  STOP

CONTAINS


  SUBROUTINE getarguments(datalabel,T,parameterlabel)

    INTENT(INOUT) datalabel,parameterlabel,T

    CHARACTER (LEN=100) :: datalabel,parameterlabel
    INTEGER :: T
    INTEGER :: counter
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
       READ(arg, '(i20)') T
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, parameterlabel)
    END IF

  END SUBROUTINE getarguments

  ! -----------------------------------------------------------------
  SUBROUTINE getsettings(datalabel, Ny)

    INTENT(IN) datalabel
    INTENT(OUT) Ny

    INTEGER :: Ny
    CHARACTER (LEN=100) :: datalabel
    CHARACTER (LEN=100) :: filename
    ! INTEGER :: dummy

    filename = trim(datalabel) // '.settings.txt'
    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

    READ(4,'(t5,i5)') Ny
    ! READ(4,'(t5,i5)') dummy ! T
    ! READ(4,'(t11,i5)') firstINT
    CLOSE(UNIT=4)

  END SUBROUTINE getsettings

  ! -----------------------------------------------------------------
  SUBROUTINE readdata(y,filename,Ny,T)
    IMPLICIT NONE

    INTENT(IN) :: filename,Ny,T
    INTENT(INOUT) :: y
    CHARACTER (LEN=200) :: filename
    CHARACTER (LEN=500) :: fmtstr

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


END PROGRAM main
! -----------------------------------------------------------------


! @\newpage\subsection{particlefilter}@
SUBROUTINE particlefilter(T, Ny, y, yNaN, Nparticles, XBAR, PARTICLEweights, DRAWllf, DRAWlike, DRAWxhat, DRAWxsig, DRAWxgain, Nstates, Nx, Nw, Ex0, sqrtVx00, DRAWsvol, Nsv, Ehbar0, Vhbar0, hbarInno, hgap0, hgaprho, sqrtVhgap, f, p, shockslopes, Nshockslopes, VSLstream,timer)

  ! use embox
  use gibbsbox, only : drawNDXsysresample
  use statespacebox, only : DLYAP
  use blaspack, only : pi, vech, ivech, eye

  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  INTENT(INOUT) :: VSLstream, timer
  INTENT(OUT)   :: XBAR, PARTICLEweights, DRAWllf, DRAWlike, DRAWxhat, DRAWxsig, DRAWxgain, DRAWsvol
  INTENT(IN)    :: T,Ny,y,yNaN, Nstates, Nx, Nw, Ex0, sqrtVx00, Nsv, Nparticles, Ehbar0, Vhbar0, hbarInno, hgap0, hgaprho, sqrtVhgap, f, p, shockslopes, Nshockslopes

  INTEGER :: J, I, K, T, Nparticles, Nstates, Nx, Ny, Nsv, p, Nw, Nshockslopes, Nsigmax

  ! OPEN MP
  ! INTEGER :: NTHREADS !, TID

  
  type(progresstimer) :: timer

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(0:T,Nparticles) :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(T,Nparticles) :: DRAWllf, DRAWlike
  DOUBLE PRECISION, DIMENSION(Nstates,T,Nparticles) :: DRAWxhat, DRAWxsig
  DOUBLE PRECISION, DIMENSION(Nstates,Ny,T,Nparticles) :: DRAWxgain
  DOUBLE PRECISION, DIMENSION(Nsv,0:T,Nparticles) :: DRAWsvol
  
  ! KBAR computation
  DOUBLE PRECISION :: XX(Nx), XBAR(Nstates,T), KKprime(Ny,Nx)

  INTEGER, DIMENSION(Nstates) :: diagndx

  ! particles
  DOUBLE PRECISION :: xposterior(Nx,Nparticles), vecSigmaX(Nx * (Nx + 1) / 2,Nparticles), zdraws(Nsv, Nparticles), hbar(Nparticles), hgap(Ny,Nparticles), SVol(Nsv,Nparticles), llf(Nparticles), likellf(Nparticles), Kprime(Ny,Nx,Nparticles), AS(Nx,Nx)
  INTEGER :: ndx(Nparticles)
  DOUBLE PRECISION :: shufflevec(Nparticles)
  ! DOUBLE PRECISION :: SHUFFLEx(Nx,Nparticles), SHUFFLEsigmax(Nx,Nx,Nparticles), SHUFFLEhbar(Nparticles), SHUFFLEhgap(Ny, Nparticles)


  ! state space objects
  DOUBLE PRECISION :: xprior(Nx), SigmaX(Nx,Nx), logdetSigmaY, ImKC(Nx,Nx)

  DOUBLE PRECISION :: f(p * Ny * Ny), shockslopes(Nshockslopes), Ehbar0, Vhbar0, sqrtVhbar0, hbarInno, hgap0(Ny), sqrtVhgap(Ny,Ny), sqrtVhgap0(Ny,Ny), Ahgap(Ny,Ny)
  DOUBLE PRECISION :: hgaprho, hgapintercept(Ny) 

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), sqrtVx00(Nx,Nx), A(Nx,Nx), B(Nx,Nw), Bsv(Nx,Nw), C(Ny,Nx,T), ytilde(Ny)
  DOUBLE PRECISION :: ygap0variance(Ny*p,Ny*p), gapshock0loadings(Ny*p,Nw)

  DOUBLE PRECISION :: minSVh(Nsv)


  ! helper for SigmaY inversion
  DOUBLE PRECISION :: thisSigmaY(Ny,Ny), thisKprime(Ny,Nx), thisC(Ny,Nx)
  INTEGER :: rndx, Nynonan, yndx(Ny)
  INTEGER :: these, offsetslopes
  ! CHARACTER (LEN=200) :: filename

  ! VSL
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
  type (vsl_stream_state) :: VSLstream
  double precision :: uniformdraws(T)

  minSVh = log(0.001d0 ** 2)
  Nsigmax = Nx * (Nx + 1) / 2

  ! prepare state space
  ! A
  A = 0.0d0
  ! unit roots for trend
  FORALL(j=1:Ny) A(j,j) = 1.0d0
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:Ny*(p-1)) A(Ny+Ny+j,Ny+j) = 1.0d0
  END IF
  ! VAR coefficients
  A(Ny+1:Ny*2,Ny+1:Nx) = transpose(reshape(f, (/ Ny * p, Ny /)))

  ! B
  B = 0.0d0
  B(1:Ny,1) = 1.0d0
  
  FORALL (j=1:Ny) B(Ny+j,1+j) = 1.0d0 
  offsetslopes = 0
  DO i = 2, Ny
      these = i-1 
      ! slopes in row i have index offsetslopes+1:offsetslopes + these
      B(Ny+i,2:1+these) = shockslopes(offsetslopes+1:offsetslopes+these)
      offsetslopes = offsetslopes + these
  END DO

  ! C
  C         = 0.0d0
  FORALL(j=1:Ny,k=1:T) C(j,j,k) = 1.0d0
  FORALL(j=1:Ny,k=1:T) C(j,Ny+j,k) = 1.0d0
  ! prepare C for missing values
  DO k=1,T
     DO j = 1, Ny
        if (yNaN(j,k)) C(j,:,k) = 0.0d0
        if (yNaN(j,k) .AND. y(j,k) /= 0.0d0 ) then
           write (*,*) 'YNAN PATTERN DOES NOT MATCH ZEROS IN Y'
        end if
     END DO
  END DO

  ! Time 0 particles
  FORALL(k=1:Nparticles) xposterior(:,k) = Ex0
  XX = Ex0

  ! draw initial values for SV
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, zdraws, 0.0d0, 1.0d0)
  ! hbar
  sqrtVhbar0 = sqrt(Vhbar0)
  FORALL (k=1:Nparticles) hbar(k) = Ehbar0 + sqrtVhbar0 * zdraws(1,k)  
  SVol(1,:) = exp(hbar * 0.5d0)
  ! hgap
  call eye(Ahgap, hgaprho)
  call DLYAP(sqrtVhgap0, Ahgap, sqrtVhgap, Ny, Ny, errcode) ! recall: sqrtVhgap is upper-triangular draw from inverse wishart, i.e. Vhgap = sqrtVhgap * sqrtVhgap'
  if (errcode /= 0) then
     write (*,*) 'DLYAP error (hGap0 variance)', errcode
     stop 1
  end if
  ! Factorize the unconditional variance
  CALL DPOTRF('L', Ny, sqrtVhgap0, Ny, errcode) ! 'L' is OK since DLAYP returns symmetric matrix
  if (errcode /= 0) then
     write (*,*) 'DPOTRF error (hgap0 variance)', errcode
     stop 1
  end if
  hgap = zdraws(2:Nsv,:)
  call dtrmm('L', 'L', 'N', 'N', Ny, Nparticles, 1.0d0, sqrtVhgap0, Ny, hgap, Ny) ! recall: sqrtVhgap0 is lower-triangular left factor, i.e. sqrtVhgap0 * sqrtVhgap'
  FORALL (i=1:Ny,k=1:Nparticles) hgap(i,k) = hgap0(i) + hgap(i,k) 
  SVol(2:Nsv,:) = exp(hgap * 0.5d0)

  vecSigmaX       = 0.0d0
  !$OMP PARALLEL DO SHARED(vecSigmaX, A, B, p, SVol, Ny, Nx, Nw, Nsv, Nparticles, sqrtVx00) PRIVATE(gapshock0loadings, sqrtVx0, ygap0variance, SigmaX, errcode) DEFAULT(NONE)
  DO k=1,Nparticles

     sqrtVx0         = sqrtVx00
     ! Fill in unconditional variance of stationary states
     ! allow for trendshockslopes, thought they are all zero here
     gapshock0loadings                       = B(Ny+1:Nx,:)
     FORALL (i=1:Nsv) gapshock0loadings(:,i) = gapshock0loadings(:,i) * SVol(i,k) 

     CALL DLYAP(ygap0variance, A(Ny+1:Nx,Ny+1:Nx), gapshock0loadings, Ny * p, Nw, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DLYAP error (ygap0variance)', errcode
        stop 1
     end if

     ! Factorize the unconditional variance
     CALL DPOTRF('L', Ny * p, ygap0variance, Ny * p, errcode) ! 'L' is OK since DLAYP returns symmetric matrix
     if (errcode /= 0) then
        write (*,*) 'DPOTRF error (ygap0variance)', errcode
        stop 1
     end if
     ! zero out the upper triangular
     FORALL (i=1:Ny*p-1) ygap0variance(i,i+1:Ny*p) = 0.0d0
     ! fill it in
     sqrtVx0(Ny+1:Nx,Ny+1:Nx) = ygap0variance

     ! store the particle's prior variance
     SigmaX = 0.0d0
     call DSYRK('U','N',Nx,Nx,1.0d0,sqrtVx0,Nx,0.0d0,SigmaX,Nx)
     call vech(vecSigmaX(:,k), SigmaX)

  END DO
  !$OMP END PARALLEL DO 

  DRAWsvol(:,0,:) = SVol
  PARTICLEweights = 1 / dble(Nparticles)
  DRAWllf         = 0.0d0
  DRAWlike        = 0.0d0

  ! uniform draws for systematic resampling
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, T, uniformdraws, 0.0d0, 1.0d0)
  
  ! prepare index for mapping diagonal element of SigmaX within vecSigmaX
  diagndx(1) = 1
  DO i=2,Nstates
     diagndx(i) = diagndx(i-1) + i
  END DO
  
  ! TID = 0
  CALL initprogressbar(timer, 15.0d0)
  DO j=1,T

     ! print *, 'Filter Time =', j

     ! 1) Draw Particles
     errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, zdraws, 0.0d0, 1.0d0)
     ! prepare hgap innovations
     ! recall: sqrtVhgap is upper-triangular right factor of inverse wishart, i.e. Vhgap = sqrtVhgap * sqrtVhgap'
     call dtrmm('L', 'U', 'N', 'N', Ny, Nparticles, 1.0d0, sqrtVhgap, Ny, zdraws(2:Nsv,:), Ny)


     ! compute Nynonan
     Nynonan = 0
     yndx = 0
     DO i=1,Ny
        if (.NOT. yNaN(i,j)) then
           Nynonan = Nynonan + 1
           yndx(Nynonan) = i
        end if
     END DO
     
     hgapintercept = (1 - hgaprho) * hgap0

     !$OMP PARALLEL DO SHARED(xposterior, vecSigmaX, SVol, hbar, hgap, Kprime, zdraws, llf, likellf, Nparticles, Nsv, hbarinno, hgapintercept, hgaprho, j, y, yNaN, A, B, C, Ny, Nx, Nw,yndx, Nynonan, minSVh) PRIVATE(ytilde, rndx, Bsv, thisSigmaY, logdetSigmaY, errcode, ImKC, xprior, SigmaX, AS, thisC, thisKprime) DEFAULT(NONE) SCHEDULE(STATIC)

     DO k = 1,Nparticles

        xprior      = 0.0d0
        SigmaX      = 0.0d0
        call ivech(SigmaX, vecSigmaX(:,k))

        hbar(k)   = hbar(k) + zdraws(1,k) * hbarinno
        SVol(1,k) = exp(hbar(k) * 0.5d0)

        hgap(:,k)     = hgapintercept + hgaprho * hgap(:,k) + zdraws(2:Nsv,k)
        SVol(2:Nsv,k) = exp(hgap(:,k) * 0.5d0)

        ! FORALL (i=1:Nsv) 
        !    WHERE (h(i,:) < minSVh(i)) h(i,:) = minSVh(i)
        ! END FORALL


        ! 2) Fill Particles into state space
        ! - SV
        FORALL (i=1:Nsv) Bsv(:,i) = B(:,i) * SVol(i,k) 

        ! 3) Kalman Filter: Prior

        ! xprior = A * xposterior(-1)
        ! call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,xprior,1)
        ! exploit trend cycle structure
        xprior(1:Ny) = xposterior(1:Ny,k)
        call DGEMV('n',Nx-Ny,Nx-Ny,1.0d0,A(Ny+1:Nx,Ny+1:Nx),Nx-Ny,xposterior(Ny+1:Nx,k),1,0.0d0,xprior(Ny+1:Nx),1)


        ! SigmaX = A * SigmaX(-1) * A' + Bsv * Bsv'
        call DSYMM('R','U', Nx, Nx,1.0d0,SigmaX,Nx,A,Nx,0.0d0,AS,Nx)
        call DGEMM('N','T',Nx,Nx,Nx,1.0d0,AS,Nx,A,Nx,0.0d0,SigmaX,NX)
        call DSYRK('U','N',Nx,Nw,1.0d0,Bsv,Nx,1.0d0,SigmaX,Nx)

       IF (Nynonan > 0) THEN

          ytilde = 0.0d0
          thisC  = 0.0d0
          FORALL (rndx = 1:Nynonan)                   ytilde(rndx)      = y(yndx(rndx),j)
          FORALL (rndx = 1:Nynonan)                   thisC(rndx,:)     = C(yndx(rndx),:,j)

          call DGEMV('n',Nynonan,Nx,-1.0d0,thisC(1:Nynonan,:),Nynonan,xprior,1,1.0d0,ytilde(1:Nynonan),1)



          ! thisKprime = thisC * SigmaX
          call DSYMM('R','U',Nynonan,Nx,1.0d0,SigmaX,Nx,thisC(1:Nynonan,:),Nynonan,0.0d0,thisKprime(1:Nynonan,:),Nynonan)
          ! thisSigmaY = thisKprime * thisC'
          call DGEMM('N','T',Nynonan,Nynonan,Nx,1.0d0,thisKprime(1:Nynonan,:),Nynonan,thisC(1:Nynonan,:),Nynonan,0.0d0,thisSigmaY(1:Nynonan,1:Nynonan),Nynonan)

          ! choleski SigmaY
          call DPOTRF('U', Nynonan, thisSigmaY(1:Nynonan,1:Nynonan), Nynonan, errcode)
          IF (errcode /= 0) THEN
             WRITE(*,'(a,i2,a)') "Error in Choleski Factorization of thisSigmaY, (DPOTRF errcode = ", errcode, ") [PARTICLEFILER]"
             print *, 'iteration j', j
             STOP 1
          END IF

          ! thisKprime = inv(SigmaY) * thisKprime
          call DPOTRS('U', Nynonan, Nx, thisSigmaY(1:Nynonan,1:Nynonan), Nynonan, thisKprime(1:Nynonan,:), Nynonan, errcode)
          IF (errcode /= 0) WRITE(*,'(a,i2,a)') "Error in Solving KPRIME (DPOTRS errcode = ", errcode, ") [PARTIFLEFILTER]"

          ! Posterior Mean and Variance
          ! xposterior = xprior + K * ytilde
          xposterior(:,k) = xprior
          call DGEMV('T',Nynonan,Nx,1.0d0,thisKprime(1:Nynonan,:),Nynonan,ytilde(1:Nynonan),1,1.0d0,xposterior(:,k),1)

          ! store thisKprime in Kprime
          Kprime(:,:,k) = 0.0d0
          FORALL (rndx = 1:Nynonan) Kprime(yndx(rndx),:,k) = thisKprime(rndx,:)


          ! ImKC = I - K * C
          call eye(ImKC)
          call DGEMM('T','N',Nx,Nx,Nynonan,-1.0d0,thisKprime(1:Nynonan,:),Nynonan,thisC(1:Nynonan,:),Nynonan,1.0d0,ImKC,Nx)
          ! SigmaX = ImKC * SigmaX * ImKC' (to preserve symmetry, using the sandwich)
          call DSYMM('R','U', Nx, Nx,1.0d0,SigmaX,Nx,ImKC,Nx,0.0d0,AS,Nx)
          call DGEMM('N','T',Nx,Nx,Nx,1.0d0,AS,Nx,ImKC,Nx,0.0d0,SigmaX,Nx)



          ! 5) Weights and LLF
          ! log(det(SigmaY))
          logdetSigmaY = 0.0d0
          DO i=1,Nynonan
             logdetSigmaY = logdetSigmaY + log(thisSigmaY(i,i))
          END DO
          logdetSigmaY = 2.0d0 * logdetSigmaY
          
          call dtrsv('U', 'T', 'N', Nynonan, thisSigmaY(1:Nynonan,1:Nynonan), Nynonan, ytilde(1:Nynonan), 1)

          ! llf
          llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde(1:Nynonan) ** 2))

          ! let's rely on the first variables always being around as long as Nynonan is true... 
          ! NOTE: thisSIGMA is transpose of left Choleski factor, so thisSigmaY(1,1) works just fine
          likellf(k)   = -0.5d0 * (log(2.0d0 * pi) + log(thisSigmaY(1,1)) + ytilde(1) ** 2)

          ! weights
          ! wstar(k) = exp(llf(k))

       ELSE
          
          xposterior(:,k)        = xprior
          ! SigmaXposterior        = SigmaXprior

          Kprime(:,:,k)          = 0.0d0
          llf(k)                 = 0.0d0
          likellf(k)             = 0.0d0

       END IF ! Nynonan > 0

       call vech(vecSigmaX(:,k), SigmaX)

     END DO ! k particles
     !$OMP END PARALLEL DO 


     ! Store NON-reweighted statistics
     DRAWllf(j,:)      = llf
     DRAWlike(j,:)     = likellf
     DRAWsvol(:,j,:)   = SVol

     FORALL(i=1:Nstates) DRAWxhat(i,j,:)      = xposterior(i,:)
     FORALL(i=1:Nstates) DRAWxsig(i,j,:)      = sqrt(vecSigmaX(diagndx(i),:))
     FORALL(i=1:Nstates) DRAWxgain(i,:,j,:)   = Kprime(:,i,:)



     if (Nynonan > 0) then

        ! Reweight particles for next round   
        PARTICLEweights(j,:) = exp(llf)
        PARTICLEweights(j,:) = PARTICLEweights(j,:) / sum(PARTICLEweights(j,:))

        call drawNDXsysresample(ndx, Nparticles, PARTICLEweights(j,:), Nparticles, uniformdraws(j))


        DO i=1,Nx
           FORALL(k=1:Nparticles) shufflevec(k) = xposterior(i,ndx(k))
           xposterior(i,:) = shufflevec
        END DO

        DO i=1,Nsigmax
              FORALL(k=1:Nparticles) shufflevec(k) = vecSigmaX(i,ndx(k))
              vecSigmaX(i,:) = shufflevec
        END DO
        

        FORALL(k=1:Nparticles) shufflevec(k) = hbar(ndx(k))
        hbar = shufflevec

        DO i=1,Ny
           FORALL(k=1:Nparticles) shufflevec(k) = hgap(i,ndx(k))
           hgap(i,:) = shufflevec
        END DO

     else ! i.e. Nynonan == 0
        PARTICLEweights(j,:) = PARTICLEweights(j-1,:) 
     end if

     
     ! compute XX
     ! 1: update XX
     xprior(1:Ny) = XX(1:Ny)
     call DGEMV('n',Nx-Ny,Nx-Ny,1.0d0,A(Ny+1:Nx,Ny+1:Nx),Nx-Ny,XX(Ny+1:Nx),1,0.0d0,xprior(Ny+1:Nx),1)
     ! 2: ytilde
     ytilde = y(:,j)
     call DGEMV('n',Ny,Nx,-1.0d0,C(:,:,j),Ny,xprior,1,1.0d0,ytilde,1)
     ! 3: integrate over gain
     forall (i=1:Ny,k=1:Nx) KKprime(i,k) = sum(PARTICLEweights(j,:) * Kprime(i,k,:))
     ! 4: xbar = xprior + K ytilde

     XX = xprior
     call DGEMV('T',Ny,Nx,1.0d0,KKprime,Ny,ytilde,1,1.0d0,XX,1)

     FORALL(i=1:Nstates) XBAR(i,j) = XX(i)

     CALL progressbarcomment(dble(j) / dble(T), timer, 'Particle Step')

  END DO ! j=1,T

END SUBROUTINE particlefilter

! @\newpage\subsection{particlesmoother}@
SUBROUTINE particlesmoother(T, Nsmoother, Nsv, SMOOTHERsvol, Nparticles, PARTICLEweights, FILTERsvol, hbarinno, Ngap, hgap0, hgaprho, sqrtVhgap, VSLstream, timer)

  ! use embox
  use gibbsbox, only: drawNDXsysresample
  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  INTENT(INOUT) :: VSLstream, timer
  INTENT(OUT)   :: SMOOTHERsvol
  INTENT(IN)    :: FILTERsvol, PARTICLEweights, T, Nsmoother, Nparticles, Nsv, hbarinno, Ngap, hgap0, hgaprho, sqrtVhgap

  INTEGER :: T, Nparticles, Nsmoother, Nsv, Ngap
  DOUBLE PRECISION :: hbarinno, hgap0(Ngap), hgaprho, hgapintercept(Ngap), sqrtVhgap(Ngap,Ngap)
  type (vsl_stream_state) :: VSLstream
  type(progresstimer) :: timer

  INTEGER :: k,i,j,errcode

  DOUBLE PRECISION, DIMENSION(0:T,Nparticles)     :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(Nsv,0:T,Nparticles) :: FILTERsvol
  DOUBLE PRECISION, DIMENSION(Nsv,0:T,Nsmoother)  :: SMOOTHERsvol

  ! filter particles
  DOUBLE PRECISION :: w(Nparticles), h(Nsv,Nparticles), hgapshock(Ngap,Nparticles)

  ! smoother objects
  DOUBLE PRECISION :: udraw(Nsmoother)
  INTEGER :: ndx(Nsmoother)
  DOUBLE PRECISION :: cdf

  INTEGER, PARAMETER :: VSLmethodUniform = 0

  ! CHARACTER (LEN=200) :: filename

  ! t=T
  ! call drawNDXpdf(ndx, Nsmoother, PARTICLEweights(T,:), Nparticles, VSLstream)
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, 1, udraw(1), 0.0d0, 1.0d0)
  call drawNDXsysresample(ndx, Nsmoother, PARTICLEweights(T,:), Nparticles, udraw(1))

  h = 2.0d0 * log(FILTERsvol(:,T,:))

  FORALL(k=1:Nsmoother) SMOOTHERsvol(:,T,k) = h(:,ndx(k))

  hgapintercept = (1.0d0 - hgaprho) * hgap0

  CALL initprogressbar(timer, 15.0d0)
  DO j=T-1,0,-1

     h = 2.0d0 * log(FILTERsvol(:,j,:))

     ! draw uniforms
     errcode = vdrnguniform(VSLmethodUniform, VSLstream, Nsmoother, udraw, 0.0d0, 1.0d0)

 
     !$OMP PARALLEL DO SHARED(PARTICLEweights, SMOOTHERsvol, h, hbarinno, udraw, hgapintercept, sqrtVhgap) PRIVATE(w,cdf,i,hgapshock) 
     DO k=1,Nsmoother

        ! draw smoothed particle k
        ! a) compute kernel weights
        ! pdf for hbar
        w = - 0.5d0 * ((SMOOTHERsvol(1,j+1,k) - h(1,:)) / hbarinno ) ** 2 

        FORALL (i=1:Nparticles) hgapshock(:,i) = SMOOTHERsvol(2:Nsv,j+1,k) - hgaprho * h(2:Nsv,i) - hgapintercept

        ! recall: sqrtVhgap is upper-triangular draw from inverse wishart, i.e. Vhgap = sqrtVhgap * sqrtVhgap'
        call  DTRSM('L', 'U', 'N', 'N', Ngap, Nparticles, 1.0d0, sqrtVhgap, Ngap, hgapshock, Ngap)
        w = w - 0.5d0 * sum(hgapshock ** 2, 1) 

        w = PARTICLEweights(j,:) * exp(w) ! works also for j=0 since PARTICLEweights = 1/ Nparticles in that case 

        w = w / sum(w)

        ! b) draw using kernel weights
        i   = 0
        cdf = 0.0d0
        DO WHILE (cdf < udraw(k))
           i   = i + 1
           cdf = cdf + w(i)
        END DO
        SMOOTHERsvol(:,j,k) = h(:,i)
        
     
     END DO ! k
     !$OMP END PARALLEL DO 

     CALL progressbarcomment(dble(T-j) / dble(T), timer, 'Smoothing Particles')
     
  END DO ! j=1,T


  ! convert log-variances back into vols
  SMOOTHERsvol = exp(SMOOTHERsvol * 0.5d0)

END SUBROUTINE particlesmoother

! @\newpage\subsection{particleSmootherX}@
SUBROUTINE particleSmootherX(T, Ny, y, yNaN, Nsv, Nsmoother, SMOOTHERsvol, NsmootherX, SMOOTHERx, Nstates, Nx, Nw, Ex0, sqrtVx00, f, p, shockslopes, Nshockslopes)

  ! use embox
  ! use gibbsbox
  use statespacebox, only : DLYAP, samplerA3B3C3nanscalar

  use vslbox
  use timerbox
  use omp_lib

  IMPLICIT NONE

  INTENT(INOUT) :: SMOOTHERx
  INTENT(IN)    :: T, Ny, Ex0, sqrtVx00, y, Nsv, Nsmoother, SMOOTHERsvol, NsmootherX, Nstates, Nx, Nw, f, p, shockslopes, Nshockslopes

  INTEGER :: T, Ny, Nsv, Nsmoother, NsmootherX, Nstates, Nx, Nw, p, Nshockslopes

  INTEGER :: j,k,i

  type(progresstimer) :: timer

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nsv,0:T,Nsmoother) :: SMOOTHERsvol
  DOUBLE PRECISION, DIMENSION(Nstates,0:T,NsmootherX,Nsmoother) :: SMOOTHERx

  ! particles
  DOUBLE PRECISION :: x(Nx,0:T), xshock(Nx,T)
  
  ! state space objects
  DOUBLE PRECISION :: f(p * Ny * Ny), shockslopes(Nshockslopes)
  DOUBLE PRECISION :: Ex0(Nx), sqrtVx00(Nx,Nx), sqrtVx0(Nx,Nx), A(Nx,Nx,T), B(Nx,Nw), Bsv(Nx,Nw,T), C(Ny,Nx,T)
  DOUBLE PRECISION :: ygap0variance(Ny*p,Ny*p), gapshock0loadings(Ny*p,Nw),dummy(Nx,Nx)

  INTEGER :: these, offsetslopes

  ! VSL
  TYPE (VSL_STREAM_STATE) :: VSLstream
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OMP
  INTEGER :: TID, NTHREADS


  ! prepare state space
  ! A
  A = 0.0d0
  ! unit roots for trend
  FORALL(j=1:Ny,k=1:T) A(j,j,k) = 1.0d0
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:Ny*(p-1),k=1:T) A(Ny+Ny+j,Ny+j,k) = 1.0d0
  END IF
  ! VAR coefficients
  FORALL(k=1:T)  A(Ny+1:Ny*2,Ny+1:Nx,k) = transpose(reshape(f, (/ Ny * p, Ny /)))

  ! B
  B = 0.0d0
  B(1:Ny,1) = 1.0d0
  
  FORALL (j=1:Ny) B(Ny+j,1+j) = 1.0d0 
  offsetslopes = 0
  DO i = 2, Ny
      these = i-1 
      ! slopes in row i have index offsetslopes+1:offsetslopes + these
      B(Ny+i,2:1+these) = shockslopes(offsetslopes+1:offsetslopes+these)
      offsetslopes = offsetslopes + these
  END DO

  ! C
  C         = 0.0d0
  FORALL(j=1:Ny,k=1:T) C(j,j,k) = 1.0d0
  FORALL(j=1:Ny,k=1:T) C(j,Ny+j,k) = 1.0d0
  ! prepare C for missing values
  DO k=1,T
     DO j = 1, Ny
        if (yNaN(j,k)) C(j,:,k) = 0.0d0
        if (yNaN(j,k) .AND. y(j,k) /= 0.0d0 ) then
           write (*,*) 'YNAN PATTERN DOES NOT MATCH ZEROS IN Y'
        end if
     END DO
  END DO

  sqrtVx0 = sqrtVx00 

  !$OMP PARALLEL SHARED(y, yNaN, SMOOTHERx, A, B, C, Ex0, Ny, p, Nx, Nw, T, SMOOTHERsvol, Nsv, Nstates, NsmootherX, Nsmoother) PRIVATE(x,Bsv,gapshock0loadings,ygap0variance,xshock,errcode,TID,NTHREADS,timer,VSLstream,dummy,i,j) FIRSTPRIVATE(sqrtVx0) DEFAULT(NONE)

  NTHREADS = 1
  TID = 0
  !$ TID = OMP_GET_THREAD_NUM()
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  errcode = vslnewstream(VSLstream, vsl_brng_mt2203 + TID + 1, 0)  
  if (errcode /= 0) then
     print *,'VSL new stream failed'
     stop 1
  end if
  ! WRITE(*,'(a16, i2, a5, i2, a8, i20, i20)') 'LAUNCHING WORKER', TID,  ' Stream ', VSLstream%descriptor1, VSLstream%descriptor2
  
  if (TID == 0) call initprogressbar(timer, 15.0d0)

  !$OMP DO 
  DO k=1,Nsmoother

     ! if (TID == 0) print *, 'k', k

     ! Factorize the unconditional variance

     ! 1) SQRTVX0
     ! a) Fill in unconditional variance of stationary states
     ygap0variance             = 0.0d0
     gapshock0loadings         = B(Ny+1:Nx,:)
     FORALL (i=1:Nsv) gapshock0loadings(:,i) = B(Ny+1:Nx,i) * SMOOTHERsvol(i,0,k)

     CALL DLYAP(ygap0variance, A(Ny+1:Nx,Ny+1:Nx,1), gapshock0loadings, Ny * p, Nw, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DLYAP error (ygap0variance)', errcode
        stop 1
     end if
     ! b) Factorize the unconditional variance
     CALL DPOTRF('L', Ny * p, ygap0variance, Ny * p, errcode) ! 'L' is OK since DLYAP returns dense-symmetrix matrix
     if (errcode /= 0) then
        write (*,*) 'DPOTRF error (ygap0variance)', errcode
        stop 1
     end if
     ! zero out the upper triangular
     FORALL (i=1:Ny*p-1) ygap0variance(i,i+1:Ny*p) = 0.0d0
     ! fill it in
     sqrtVx0(Ny+1:Nx,Ny+1:Nx) = ygap0variance

     ! 2) BSV 
     FORALL (j=1:T,i=1:Nsv) Bsv(:,i,j)   = B(:,i) * SMOOTHERsvol(i,j,k)
     
     ! 3) KALMAN SMOOTHER: NsmootherX draws
     DO j=1,NsmootherX
        ! call samplerA3B3C3charly(x(:,j,:),xshock,y,T,Ny,Nx,Nw,A,Bsv,C,Ex0,sqrtVx0,VSLstreams(TID),errcode)
        ! note: xshock is a dummy

        call samplerA3B3C3nanscalar(x,xshock,dummy,y,yNaN,T,Ny,Nx,Nw,A,Bsv,C,Ex0,sqrtVx0,VSLstream,errcode)
        FORALL (i=1:Nstates) SMOOTHERx(i,:,j,k)  = x(i,:)

     END DO ! j


     ! timer assumes that equal chunks are spread across workers
     if (TID == 0) call progressbarcomment(dble(k * NTHREADS) / dble(Nsmoother), timer, 'Smoothing Particles X')
  END DO ! k
  !$OMP END DO 


  errcode = vsldeletestream(VSLstream)     


  !$OMP END PARALLEL


END SUBROUTINE particleSmootherX
