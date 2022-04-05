PROGRAM main

  USE embox, only : hrulefill, savemat, savevec, storeEstimates, storeEstimatesOMP, loft, timestampstr, es30d16, int2str, mean
  USE blaspack, only : eye, vech
  USE gibbsbox, only : igammaDraw, iwishDraw, GelmanTest1, simpriormaxroot, ineffparzen, ineffbrtltt, ineffbatch
  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  INTEGER :: Ny, p, Nbar, Ngap, Nstates, Nx, Nsv, Nf, Nsigmagap
  
  LOGICAL, PARAMETER :: doTimestamp = .false.

  ! double precision :: sqrtVf0_firstlag, sqrtVf0_general
  DOUBLE PRECISION :: minnesotaTheta1Vf, minnesotaTheta2Vf, minnesotaTheta3Vf

  INTEGER :: Nsim, Burnin, Nstreams
  INTEGER :: T,j,k,i,status,Ndraws

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWsvol, DRAWstates
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWf, DRAWsigmagap, DRAWmaxlambda, DRAWhvarbar

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: minSV, maxSV, Ef0
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sqrtVf0, SigmaGapT
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SigmaY
  
  DOUBLE PRECISION :: hvarbarT, Ehbar0, Vhbar0
  INTEGER :: SigmaGapDof, hvarbarDof

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: GELMANstates, GELMANsvol
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GELMANf, GELMANsigmagap, GELMANhvarbar

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: theta
  !  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: thetadraws

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ydata
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN


  TYPE(progresstimer) :: timer
  CHARACTER (LEN=200) :: filename, datafile, nandatafile, filext, datalabel

  ! VSL Random Stuff
  type (vsl_stream_state) :: VSLdefaultstream, VSLstream
  integer :: seed
  integer :: errcode

  ! OPEN MP
  INTEGER :: NTHREADS, TID

  ! ----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  ! ----------------------------------------------------------------------------------
  ! runtime parameters :start:
  ! first: set default values
  Nsim    = 10 ** 4
  Burnin  = 10 ** 5

  ! Nsim    = 10 ** 3
  ! Burnin  = 10 ** 3

  datalabel = 'INFTRM'

  ! original settings
  ! sqrtVf0_firstlag = 0.10d0
  ! sqrtVf0_general  = 0.01d0


  p = 12

  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  call hrulefill
  print *, 'MCMC estimation of PaddingtonGAPconst'
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
  call getarguments(datalabel,p,T,Nsim,burnin)
  call getsettings(datalabel,Ny)

  call hrulefill
  print *, 'datalabel: ', datalabel
  print *, 'with Minnesota VAR prior'
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

  filext = '.' // trim(datalabel) // '.T' // trim(int2str(T)) // '.gapCONST.dat'
  IF (doTimeStamp) filext = '.' // timestampstr() // filext
  filext = '.' // 'notrendslopes' // filext


  Nbar = Ny
  Ngap = Ny
  Nstates = Nbar + Ngap
  Nx = Nbar + p * Ngap
  Nsv = 1
  Nf  = Ngap * Ngap * p
  Nsigmagap = Ngap * (Ngap + 1) / 2

  datafile    = trim(datalabel) // '.yData.txt'
  nandatafile = trim(datalabel) // '.yNaN.txt'

  ALLOCATE (yNaN(Ny,T), ydata(Ny,T), SigmaY(Ny,Ny), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Ydata)'
  END IF

  CALL readdata(ydata,datafile,Ny,T)
  CALL readnandata(yNaN,nandatafile,Ny,T)


  ! Model parameters and priors
  ALLOCATE (minSV(Nsv), maxSV(Nsv), Ef0(Nf), sqrtVf0(Nf,Nf), SigmaGapT(Ngap,Ngap))

  print *, 'DONE ALLOCATION'


  ! SV
  minSV = 1.0d-3
  maxSV = 40.0d0 

  ! TREND SV
  Vhbar0     = 1.0d0
  Ehbar0     = log(1 / 12.0d0) - Vhbar0 * 0.25d0;
  hvarbarDof = 2 + 24
  hvarbarT   = (0.2d0 ** 2) / 3.0d0 * dble(hvarbarDof - 2)

  ! Gap VCV SigmaGap
  call eye(SigmaGapT, 0.5d0 ** 2) 
  SigmaGapDof = Ny + 2



  ! VAR coefficients
  Ef0  = 0.0d0
  ! call eye(sqrtVf0, sqrtVf0_general)
  ! ! widen the prior for first-own-lag coefficients
  ! DO j=1,NGap
  !    k = (j-1)*(Ngap * p) + j
  !    sqrtVf0(k,k) = sqrtVf0_firstlag
  ! END DO

  minnesotaTheta1Vf = 0.1d0   ! overall shrinkage ! set to 0.05 for tightPrior
  minnesotaTheta2Vf = 0.5d0   ! cross shrinkage
  minnesotaTheta3Vf = 2.0d0   ! decay

  call minnesotaccmVCVsqrt(sqrtVf0, Ny, p, minnesotaTheta1Vf, minnesotaTheta2Vf, minnesotaTheta3Vf)


  ! lambda1Vf = 0.5d0
  ! lambda2Vf = 0.2d0
  ! call minnesotaVCVsqrt(sqrtVf0, Ny, p, lambda1Vf, lambda2Vf)




  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data= ' // datalabel
  print *, 'T = ', T
  print *, 'Ny= ', Ny
  print *, 'Nsim= ', Nsim
  print *, 'Burnin= ', Burnin
  print *, 'Nstreams= ', Nstreams
  print *, 'p= ', p
  CALL HRULEFILL

  ! allocate memory for draws
  ALLOCATE (DRAWhvarbar(1,Nsim,Nstreams),DRAWmaxlambda(1,Nsim,Nstreams), DRAWf(Nf,Nsim,Nstreams), DRAWsigmagap(Nsigmagap,Nsim,Nstreams), DRAWstates(Nstates,0:T,Nsim,Nstreams), DRAWsvol(Nsv,0:T,Nsim,Nstreams), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (draws)'
  END IF

  WRITE(*,*) 'STARTING SWEEPS'

  !$OMP PARALLEL DO SHARED(DRAWstates,DRAWsvol,DRAWhvarbar,DRAWf,DRAWsigmagap,DRAWmaxlambda, Nstreams,Nsim,Burnin,T,p,ydata,yNaN,Ny,Nstates,Nx,Nsv,Ehbar0,Vhbar0,hvarbarT,hvarbarDof,minSV,maxSV,Nf,Ef0,sqrtVf0,Nsigmagap,SigmaGapT,Ngap,SigmaGapDof), PRIVATE(VSLstream,TID,errcode,timer), DEFAULT(NONE) SCHEDULE(STATIC)

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
     CALL thissampler(T,p,ydata,yNaN,Ny,DRAWstates(:,:,:,j),Nstates,Nx,DRAWsvol(:,:,:,j),Nsv,Ehbar0,Vhbar0,DRAWhvarbar(:,:,j),hvarbarT,hvarbarDof,minSV,maxSV, DRAWf(:,:,j), Nf, Ef0, sqrtVf0, DRAWmaxlambda(:,:,j),DRAWsigmagap(:,:,j), Nsigmagap, SigmaGapT, Ngap, SigmaGapDof, Nsim,Burnin,VSLstream,timer)

     ! compute and report INEFF factors 
     print *, 'TID', TID, 'f', mean(ineffparzen(DRAWf(:,:,j))), mean(ineffbrtltt(DRAWf(:,:,j))), mean(ineffbatch(DRAWf(:,:,j)))
     print *, 'TID', TID, 'a', mean(ineffparzen(DRAWhvarbar(:,:,j))), mean(ineffbrtltt(DRAWhvarbar(:,:,j))), mean(ineffbatch(DRAWhvarbar(:,:,j))) 

     print *, 'TID', TID, 'SV', mean(ineffparzen(reshape(DRAWsvol(:,:,:,j), (/Nsv * (T + 1), Nsim /)))), mean(ineffbrtltt(reshape(DRAWsvol(:,:,:,j), (/Nsv * (T + 1), Nsim /)))), mean(ineffbatch(reshape(DRAWsvol(:,:,:,j), (/Nsv * (T + 1), Nsim /)))) 
     print *, 'TID', TID, 'States', mean(ineffparzen(reshape(DRAWstates(:,:,:,j), (/Nstates * (T + 1), Nsim /)))), mean(ineffbrtltt(reshape(DRAWstates(:,:,:,j), (/Nstates * (T + 1), Nsim /)))), mean(ineffbatch(reshape(DRAWstates(:,:,:,j), (/Nstates * (T + 1), Nsim /)))) 



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
  WRITE(4,'(a20,I20)') 'Sims: ', Nsim
  WRITE(4,'(a20,I20)') 'Burnin: ', Burnin
  WRITE(4,'(a20,I20)') 'Streams: ', Nstreams
  WRITE(4,'(a20,I20)')   'p               = ', p
  WRITE(4,'(a60)') repeat('-',60)
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


  ! differences in TAU0
  DO i=2,Ny
     filename  = 'TAU0DELTA' // trim(int2str(i)) // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(1,(j-1) * Nsim + k) = DRAWstates(i,0,k,j) - DRAWstates(1,0,k,j)
     CALL savevec(theta(1,:),filename)
     WRITE (*,*) 'STORED TAU0DELTA', i
  END DO


  filename  = 'TAU' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(1,1:T,k,j)
  CALL storeEstimates(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED TAU'

  filename  = 'SIGTREND' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWsvol(1,1:T,k,j)
  CALL storeEstimates(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED SIGTREND'

  DEALLOCATE(theta)

  ALLOCATE(theta(Ndraws,1))
  filename  = 'HVARBAR' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta((j-1) * Nsim + k,:) = DRAWhvarbar(:,k,j)
  CALL savemat(theta, filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HVARBAR'

  ALLOCATE(theta(Nf,Ndraws))
  filename  = 'F' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWf(:,k,j)
  CALL storeEstimates(theta,Nf,Ndraws,filename)
  ! CALL savemat(theta, filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED F'

  ALLOCATE(theta(Nsigmagap,Ndraws))
  filename  = 'SIGMAGAP.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWsigmagap(:,k,j)
  CALL savemat(transpose(theta), filename)
  filename  = 'SIGMAGAP' // filext
  CALL storeEstimates(theta,Nsigmagap,Ndraws,filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED SIGMAGAP'

  ALLOCATE(theta(1,Ndraws))
  filename  = 'MAXLAMBDA' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWmaxlambda(:,k,j)
  ! CALL storeEstimates(theta,1,Ndraws,filename)
  CALL savevec(theta(1,:), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED MAXLAMBDA'


  WRITE(*,*) 'ALL ESTIMATES HAVE BEEN STORED (' // trim(adjustl(filext)) // ')'

  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  ! GELMAN
  ! ----------------------------------------------------------------------------
  
  ALLOCATE (GELMANstates(Nstates,T), GELMANsvol(Nsv,T), GELMANf(Nf), GELMANhvarbar(1), GELMANsigmagap(Nsigmagap), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Gelman statistics)'
  END IF

  CALL HRULEFILL
  WRITE (*,*) 'GELMAN STORAGE ALLOCATED'
  CALL HRULEFILL


  !$OMP PARALLEL SHARED(DRAWstates,DRAWsvol,GELMANstates,GELMANsvol,DRAWf,GELMANF,DRAWsigmagap,GELMANSIGMAGAP)

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
  DO j = 1, Nf
     call GelmanTest1(GELMANf(j), DRAWf(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO
  
  !$OMP DO 
  DO j = 1, Nsigmagap
     call GelmanTest1(GELMANsigmagap(j), DRAWsigmagap(j,:,:), Nsim, Nstreams)
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


  IF (ALL(ABS(GELMANhvarbar - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HVARBAR'
  ELSE
     WRITE(*,*) 'HVARBAR: GELMAN FAILURE, Max SRstat=', maxval(GELMANhvarbar)
  END IF
  filename = 'GELMAN.hvarbar' // filext
  call savevec(GELMANhvarbar, filename)


  IF (ALL(ABS(GELMANf - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for F'
  ELSE
     WRITE(*,*) 'F: GELMAN FAILURE, Max SRstat=', maxval(GELMANf)
  END IF
  filename = 'GELMAN.f' // filext
  call savevec(GELMANf, filename)

  IF (ALL(ABS(GELMANsigmagap - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for SIGMAGAP'
  ELSE
     WRITE(*,*) 'SIGMAGAP: GELMAN FAILURE, Max SRstat=', maxval(GELMANsigmagap)
  END IF
  filename = 'GELMAN.sigmagap' // filext
  call savevec(GELMANsigmagap, filename)

  DEALLOCATE (GELMANstates, GELMANsvol, GELMANf, GELMANsigmagap)


  ! ----------------------------------------------------------------------------
  ! FINISHED: GELMAN
  ! ----------------------------------------------------------------------------

  call hrulefill
  WRITE(*,*) 'DONE, BYE BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  ! ----------------------------------------------------------------------------
  ! CLEANUP
  ! ----------------------------------------------------------------------------

  DEALLOCATE (DRAWmaxlambda, DRAWhvarbar, DRAWstates,DRAWsvol, DRAWf, DRAWsigmagap, STAT=status)
  DEALLOCATE (ydata, yNaN)
  DEALLOCATE (minSV, maxSV, Ef0, sqrtVf0, SigmaGapT)
  ! VSLstreams
  errcode = vsldeletestream(VSLdefaultstream)     
  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill
  STOP

CONTAINS


  SUBROUTINE minnesotaccmVCVsqrt(sqrtVf0, N, p, theta1, theta2, theta3)

    integer, intent(in) :: N, p
    double precision, intent(inout), dimension(N*(N*p),N*(N*p)) :: sqrtVf0
    double precision, intent(in) :: theta1, theta2, theta3

    integer :: ndxVec, ndxLHS, thislag, ndxRHS
    integer :: Nf

    ! init
    Nf = N * N * p
    call eye(sqrtVf0) !  = 0.0d0

    ! construct prior variance
    ndxVec = 0
    do ndxLHS = 1,N

       ! lags
       do thislag = 1,p
          do ndxRHS = 1,N
             ndxVec = ndxVec + 1
             if (ndxLHS == ndxRHS) then
                ! own lag
                sqrtVf0(ndxVec,ndxVec) = sqrt(theta1 / (dble(thislag) ** theta3))
             else
                ! cross lags
                sqrtVf0(ndxVec,ndxVec) = sqrt(theta1 * theta2 /  (dble(thislag) ** theta3))
             end if
          end do
       end do
    end do

  END SUBROUTINE minnesotaccmVCVsqrt

  SUBROUTINE minnesotaVCVsqrt(sqrtVf0, N, p, lambda1, lambda2)

    integer, intent(in) :: N, p
    double precision, intent(inout), dimension(N*N*p,N*N*p) :: sqrtVf0
    double precision, intent(in) :: lambda1, lambda2

    integer :: ndxVec, ndxLHS = 1, ndxLag = 1, ndxRHS = 1
    integer :: Nf

    Nf = N * N * p
    sqrtVf0 = 0.0d0

    ndxVec = 0
    do ndxLHS = 1,N
       do ndxLag = 1,p
          do ndxRHS = 1,N

             ndxVec = ndxVec + 1

             if (ndxLHS == ndxRHS) then
                sqrtVf0(ndxVec,ndxVec) = lambda1 / dble(ndxLag)
             else
                sqrtVf0(ndxVec,ndxVec) = lambda1 * lambda2 / dble(ndxLag)
             end if
             
          end do
       end do
    end do
    
  END SUBROUTINE minnesotaVCVsqrt

  ! -----------------------------------------------------------------
  SUBROUTINE getarguments(datalabel,p,T,Nsim,burnin)

    INTENT(INOUT) datalabel, p, T, Nsim, burnin

    INTEGER :: T, p, Nsim, burnin
    INTEGER :: counter
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

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') Nsim
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') burnin
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

END PROGRAM main
! -----------------------------------------------------------------

SUBROUTINE thissampler(T,p,y,yNaN,Ny,DRAWstates,Nstates,Nx,DRAWsvol,Nsv,Ehbar0,Vhbar0,DRAWhvarbar,hvarbarT,hvarbarDof,minSV,maxSV, DRAWf, Nf, Ef0, sqrtVf0, DRAWmaxlambda,DRAWsigmagap, Nsigmagap, SigmaGapT, Ngap, SigmaGapDof, Nsim,Burnin,VSLstream,timer)

  use embox, only : hrulefill, savemat, savevec
  use blaspack, only : eye, vech
  use gibbsbox, only : varmaxroot, igammadraw, drawRW, iwishcholdraw, bayesVARbarshock, bayesVAR, stochvolKSC0, varianceDRAW, vcvcholDrawTR
  use statespacebox, only : DLYAP, samplerA3B3C3nanscalar
  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  INTENT(INOUT) :: DRAWstates,DRAWsvol,DRAWhvarbar,DRAWf,DRAWsigmagap, VSLstream, timer
  INTENT(IN)    :: T,p,y,yNaN,Ny,Nstates,Nx,Nsv,Ehbar0,Vhbar0,hvarbarT,hvarbarDof,minSV,maxSV,Nsim,Burnin, Nf, Ef0, sqrtVf0,Nsigmagap, SigmaGapT, SigmaGapDof, Ngap

  INTEGER :: J, I, T, Nsim, Burnin, TotalSim, thisdraw, Nstates, Nx, Ny, Nsv, p, Nbar, Ngap, Nf, status, SigmaGapDof, Nsigmagap, Nw

  type (vsl_stream_state) :: VSLstream
  type(progresstimer) :: timer

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nstates,0:T,Nsim) :: DRAWstates
  DOUBLE PRECISION, DIMENSION(Nsv,0:T,Nsim) :: DRAWsvol
  DOUBLE PRECISION, DIMENSION(Nf,Nsim) :: DRAWf
  DOUBLE PRECISION, DIMENSION(1,Nsim) :: DRAWmaxlambda, DRAWhvarbar
  DOUBLE PRECISION, DIMENSION(Nsigmagap,Nsim) :: DRAWsigmagap

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), A(Nx,Nx,T), B(Nx,Nsv + Ngap,T), C(Ny,Nx,T), Ef0(Nf), sqrtVf0(Nf,Nf), iVf0(Nf,Nf), f(Nf), SigmaGapT(Ngap,Ngap), SigmaGap(Ngap,Ngap), invSigmaGap(Ny,Ny), maxlambda, Iy(Ny,Ny),SigmaStarT(Nx,Nx)

  DOUBLE PRECISION, PARAMETER :: unity = 1.0d0 - 1.0d-4

  DOUBLE PRECISION :: gapdraw(-(p-1):T,Ngap), gapshock(1:T,Ngap), gap0variance(Ngap*p,Ngap*p)

  DOUBLE PRECISION, DIMENSION(Nsv) :: minSV, maxSV
  DOUBLE PRECISION :: Ehbar0,Vhbar0,hvarbarT, hvarbar
  INTEGER :: hvarbarDof

  DOUBLE PRECISION, DIMENSION(Nx, 0:T) :: x
  DOUBLE PRECISION, DIMENSION(Nx, T) :: xshock !, xdisturbance
  DOUBLE PRECISION, DIMENSION(T) :: TrendInno 
  DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: SVol
  DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: h
  DOUBLE PRECISION, DIMENSION(Nsv,1:T) :: hshock

  ! stack management
  integer, parameter :: stacksize = 5, maxShakes = 1
  double precision, dimension(Nsv,0:T,stacksize) :: PREV_SVOL
  double precision, dimension(Nf,stacksize) :: PREV_F
  double precision, dimension(1,stacksize) :: PREV_HVARBAR
  double precision, dimension(Ngap,Ngap,stacksize) :: PREV_SIGMAGAP


  integer, dimension(stacksize) :: stackRegister
  integer :: lastInStack, shakes, stackResetCount
  logical :: OK 
  CHARACTER (LEN=200) :: resetmsg

  ! VSL
  INTEGER :: errcode

  stackRegister = -1 ! dummy value for unregistered stack
  stackResetCount = -1

  TotalSim = Burnin + Nsim ! only used for updating the progress bar

  Nbar = Ny
  Nw = Nsv + Ngap
  call eye(Iy)

  ! prepare state space
  Ex0 = 0.0d0
  Ex0(1:Nbar) = 2.0d0
  ! call eye(sqrtVx0, 1.0d2)
  ! sqrtVx0, expressed as lower  triangular-choleski factor (sqrtVx0 * sqrtVx0')
  sqrtVx0 = 0.0d0
  sqrtVx0(1:Nbar,1) = 100.d0                 ! uncertainty about "primary" trend
  FORALL (j=2:Nbar) sqrtVx0(j,j) = 2.0d0     ! uncertainty about other trend given primary

  A = 0.0d0
  ! unit roots for bar
  FORALL (j=1:Nbar) A(j,j,:) = 1.0d0

  ! kompanion
  FORALL (j = 1 : Ngap * (p - 1))
     A(Nbar+Ngap+j,Nbar+j,:) = 1.0d0
  END FORALL

  B = 0.0d0

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


  shakes = 1
  j = 1
  lastInStack = 1 ! redundant, but avoids compiler warning
  resetmsg = 'DRAW 0'


  DO 


     ! WRITE (*,*) 'STEP ', j

     IF (j == 1) THEN

        call initprogressbar(timer,timer%wait,timer%rank)
        stackResetCount = stackResetCount + 1

        call hrulefill
        print *, 'RE-INIT STREAM ', timer%rank, '(reset count = ', stackResetCount, ', ',  trim(resetmsg), ')'
        call hrulefill

        ! init SV
        ! WRITE (*,*) 'INI SV ...'
           DO ! repeat until SV inside bounds

              ! hvar
              call igammaDraw(hvarbar, hvarbarT, hvarbarDof, VSLstream)

              CALL drawRW(h(1,:),T,sqrt(hvarbar),Ehbar0,Vhbar0, VSLstream)
              SVol(1,:) = exp(h(1,:) * 0.5d0)

              OK = ALL(SVol(1,:) >= minSV(1)) .AND. ALL(SVol(1,:) <= maxSV(1))
              IF (OK) EXIT

           END DO
        ! WRITE (*,*) 'DONE INIT SV ...'

        ! init VAR coefficients f
        DO
           errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nf, f, 0.0d0, 1.0d0)
           CALL DTRMV('U', 'T', 'N', Nf, sqrtVf0, Nf, f, 1)
           f = f + Ef0
           call VARmaxroot(maxlambda, f, Ngap, p)
           OK = maxlambda < unity
           IF (OK) EXIT
        END DO

        ! init VCV of Gaps
        if (SigmaGapDof > 0) then
           call iwishcholDraw(SigmaGap, SigmaGapT, SigmaGapDof, Ny, VSLstream)
        else
           call iwishcholDraw(SigmaGap, Iy, Ny + 2, Ny, VSLstream)
        end if
        

        ! init stack
        lastInStack = 1

        PREV_SVOL(:,:,lastInStack)      = SVol
        PREV_HVARBAR(:,lastInStack)     = hvarbar
        PREV_SIGMAGAP(:,:,lastInStack)  = SigmaGap
        PREV_F(:,lastInStack)           = f

        stackRegister(lastInStack)   = j

     ELSE ! j > 1

        SVol      = PREV_SVOL(:,:,lastInStack)
        hvarbar   = PREV_HVARBAR(1,lastInStack) 
        f         = PREV_F(:,lastInStack) 
        SigmaGap  = PREV_SIGMAGAP(:,:,lastInStack)

     END IF ! j == 1

     ! check inputs
     if (any(SVol <= 0.0d0)) then
        WRITE(*,*) 'negative SVOL?!'
        STOP 1
     end if


     ! ---------------------------------------------------------------------------
     ! SMOOTHING SAMPLER
     ! ---------------------------------------------------------------------------

     ! Fill SV into statespace
     FORALL (i=1:Ny) B(i,1,:)        = SVol(1,1:T)

     ! Fill in chol(SigmaGap) -- it is upper triangular b/o iwishdraw
     FORALL (i=1:T) B(Nbar+1:Nbar+Ngap,2:1+Ngap,i) = SigmaGap

     ! if (B(Nbar+Ngap,2,1) .ne. 0.0d0) then
     !    call savemat(B(:,:,1), 'B.debug')
     !    print *, 'houston'
     !    stop 1
     ! end if
     

     ! Fill in VAR coefficients
     FORALL (i=1:T) A(Nbar+1:Nbar+Ngap,Nbar+1:Nbar+Ngap*p,i) = transpose(reshape(f, (/ Ngap * p, Ngap /)))

     ! Fill in unconditional variance of stationary states
     CALL DLYAP(gap0variance, A(Nbar+1:Nx,Nbar+1:Nx,1), B(Nbar+1:Nx,:,1), Ngap * p, Nw, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DLYAP error (gap0variance)', errcode
        call savemat(B(Nbar+1:Nx,:,1), 'B.debug')
        call savemat(A(Nbar+1:Nx,Nbar+1:Nx,1), 'A.debug')
        stop 1
     end if
     ! Factorize the unconditional variance
     CALL DPOTRF('U', Ngap * p, gap0variance, Ngap * p, errcode)
     if (errcode /= 0) then
        write (*,*) 'DPOTRF error (ygapvariance)', errcode
        call savemat(B(Nbar+1:Nx,:,1), 'B.debug')
        call savemat(A(Nbar+1:Nx,Nbar+1:Nx,1), 'A.debug')
        stop 1
     end if

     ! zero out the lower triangular
     FORALL (i=2:Ngap * p) gap0variance(i,1:i-1) = 0.0d0
     ! fill it in
     sqrtVx0(Nbar+1:Nx,Nbar+1:Nx) = transpose(gap0variance)


     ! call savemat(A(:,:,1), 'A.debug')
     ! call savemat(B(:,:,1), 'B.debug')
     ! call savemat(sqrtVx0, 'sqrtVx0.debug')
     ! stop 33

     ! call smoothing sampler
     ! CALL samplerA3B3C3NaN(x,xshock,y,yNaN,T,Ny,Nx,Nw,A,B,C,Ex0,sqrtVx0,VSLstream)

     CALL samplerA3B3C3nanscalar(x,xshock,SigmaStarT,y,yNaN,T,Ny,Nx,Nw,A,B,C,Ex0,sqrtVx0,VSLstream, errcode)
     

     if (any(isnan(x))) then
        write(*,*) 'NaN draws from Kalman Sampler'
        stop 1
     end if


     ! ---------------------------------------------------------------------------
     ! VAR Step
     ! ---------------------------------------------------------------------------

     ! gapbarshock
     trendinno = xshock(1,:) / SVol(1,1:T)

     ! a) construct SigmaGap and then invert it
     ! exploiting factorization in SigmaGap (could also carry forward wishart draw ...
     call DTRTRI('U', 'N', Ngap, SigmaGap, Ngap, status)
     if (status /= 0) then
        write(*,*) 'DPOTRR ERROR, INFO: ', status, ' [SIGMAGAP - PADDINGTONSAMPLER]'
        stop 1
     end if
     call dsyrk('U', 'T', Ngap, Ngap, 1.0d0, SigmaGap, Ngap, 0.0d0, invSigmaGap, Ngap)

     ! b) collect lagged gaps (note: will be zero under tight prior)
     FORALL (i=1:p) gapdraw(1-i,:) = x(Nbar+(i-1)*Ngap+1:Nbar+i*Ngap,0)
     gapdraw(1:T,:)                = transpose(x(Nbar+1:Nbar+Ngap,1:T))

     if (any(abs(gapdraw) > 1000)) then
        print *, 'weird gaps (too big)', j
        call savemat(transpose(x), 'x.dat.debug')
        call savemat(transpose(y), 'y.dat.debug')
        call savemat(transpose(SVol), 'SVol.dat.debug')
        call savemat(A(:,:,1), 'A.dat.debug')
        call savevec(f, 'f.dat.debug')
        stop 2
     end if


     ! d) draw VAR coefficients
     ! recall: SigmaGap contains inverse of SigmaGap at this point
     
     call bayesVAR(f, gapshock, p, gapdraw, Ngap, T, invSigmaGap, Ef0, iVf0, VSLstream)
     

     ! if (any(isnan(f))) then
     !    write(*,*) 'NaN draws from bayesVAR'
     !    call display(SigmaGap, 'f6.2')
     !    call savemat(gapdraw, 'debug.gapdraw.dat')
     !    stop 1
     ! end if

     ! d) check stability
     call VARmaxroot(maxlambda, f, Ngap, p)


     ! print *, 'maxroot', maxlambda
     ! stop 1

     OK = maxlambda < unity
     IF (.NOT. OK) resetmsg = 'GAP VAR'

     ! ---------------------------------------------------------------------------


     IF (OK) THEN 

        ! ---------------------------------------------------------------------------
        ! SV Step: JPR
        ! ---------------------------------------------------------------------------
        CALL stochvolKSC0(h(1,:), xshock(1,:), T, sqrt(hvarbar), Ehbar0, Vhbar0, VSLstream)
        SVol(1,:) = exp(h(1,:) * 0.5d0)
        OK = ALL(SVol(1,:) >= minSV(1)) .AND. ALL(SVol(1,:) <= maxSV(1))
        IF (.NOT. OK) resetmsg = 'SVOL'

     END IF

     IF (OK) THEN ! move forward

        ! ---------------------------------------------------------------------------
        ! hvarbar draw
        ! ---------------------------------------------------------------------------
        hshock = h(:,1:T) - h(:,0:T-1)
        call varianceDraw(hvarbar, hvarbarT, hvarbarDof, hshock(1,:), T, VSLstream)


        ! ---------------------------------------------------------------------------
        ! Gap VCV Step
        ! call vcvDraw(SigmaGap, SigmaGapT, SigmaGapDof, gapshock, T, Ngap, VSLstream)
        call vcvcholDrawTR(SigmaGap, SigmaGapT, SigmaGapDof, gapshock, T, Ngap, VSLstream)
        ! ---------------------------------------------------------------------------

        ! store Draws after Burnin
        IF (j > BURNIN) THEN
           thisdraw = j - BURNIN
           DRAWstates(1:Ny,:,thisdraw) = x(1:Nbar,:)
           DRAWstates(Ny+1:Nstates,:,thisdraw) = x(Nbar+1:Nbar+Ngap,:)

           DRAWf(:,thisdraw) = f
           DRAWhvarbar(:,thisdraw) = hvarbar
           DRAWmaxlambda(:,thisdraw) = maxlambda
           DRAWsvol(:,:,thisdraw) = SVol

           call vech(DRAWsigmagap(:,thisdraw), SigmaGap)

        END IF ! J > BURNIN

        ! update stack
        lastInStack = minloc(stackRegister,1) ! find a free element on the stack
        stackRegister(lastInStack) = j + 1 ! notice the +1
        
        PREV_SVOL(:,:,lastInStack) = SVol
        PREV_HVARBAR(:,lastInStack) = hvarbar
        PREV_F(:,lastInStack) = f
        PREV_SIGMAGAP(:,:,lastInStack) = SigmaGap

        j = j + 1

     ELSE ! NOT(OK): shake or move back

        OK = .TRUE.
        IF (shakes < maxshakes) THEN
           shakes = shakes + 1
           WRITE (*,'("Stream ", i2, ": shaking at step ", i5)') timer%rank, j

        ELSE ! move back

           shakes = 1
           IF (j > 1) THEN
              ! WRITE (*,'("Stream ", i2, ": moving back at step ", i5)') timer%rank, j
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
     
           END IF ! j > 1

        END IF ! shakes 

     END IF ! OK

     if (j > TotalSim) EXIT

     CALL progressbar(dble(j) / dble(TotalSim), timer)


  END DO ! while j < TotalSim

END SUBROUTINE thissampler
