MODULE embox

  USE vslbox

  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793115997963468544d0, logtwopi = log(2.0d0 * pi), sqrttwopi = sqrt(2.0d0 * pi), invsqrttwopi = 1 / sqrt(2.0d0 * pi)


CONTAINS


! @\newpage\subsection{loft}@
  FUNCTION loft(filename)
    INTENT(IN) :: filename


    character(len=200) :: filename, dummy
    integer :: loft,status 

    open(unit=4,file=filename)
    loft = 0
    do
       read(4,*,IOSTAT=status) dummy
       IF (status /=0) EXIT
       loft = loft + 1
    end do

    close(unit=4)

  END FUNCTION loft

! @\newpage\subsection{mean}@
  PURE FUNCTION mean(x) result(m)
    INTENT(IN) :: x

    double precision, dimension(:) :: x
    double precision :: m

    m = sum(x) / size(x)

  END FUNCTION mean

! @\newpage\subsection{variance}@
  PURE FUNCTION variance(x) result(v)
    INTENT(IN) :: x

    double precision, dimension(:) :: x
    double precision :: v
    integer :: N

    N = size(x)
    v = sum((x - (sum(x) / N)) ** 2) / N

  END FUNCTION variance

! @\newpage\subsection{normpdf}@
  ELEMENTAL FUNCTION normpdf(x) 

    INTENT(IN) :: x

    double precision :: x, normpdf

    normpdf = exp(-0.5d0 * (logtwopi + x ** 2))

  END FUNCTION normpdf


! @\newpage\subsection{timestampstr}@
  FUNCTION timestampstr()
    character(18) :: timestampstr
    character(8)  :: date
    character(10) :: time
    call date_and_time(date,time)

    timestampstr = date//'d'//time(1:4)//'h'//time(5:6)//time(8:9)
  END FUNCTION timestampstr

! @\newpage\subsection{timestamp}@
  FUNCTION timestamp()
    INTEGER :: timestamp(8)
    call date_and_time(VALUES=timestamp)

  END FUNCTION timestamp

! @\newpage\subsection{HRULEFILL}@
  SUBROUTINE HRULEFILL

    WRITE (*,*) '---------------------------------------------------------------------------'

  END SUBROUTINE HRULEFILL

! @\newpage\subsection{display}@
  SUBROUTINE display(x,fs)
    IMPLICIT NONE

    INTENT(IN) :: x

    INTEGER rows,cols,j
    CHARACTER (LEN=*) :: fs
    DOUBLE PRECISION, DIMENSION(:,:) :: x

    cols = size(x,2)
    rows = size(x,1)

    WRITE(*,fmtstr(fs,cols)) (x(j,:), j=1,rows)

  END SUBROUTINE display

! @\newpage\subsection{displayvec}@
  SUBROUTINE displayvec(x,fs)
    IMPLICIT NONE

    INTENT(IN) :: x

    INTEGER rows,j
    CHARACTER (LEN=*) :: fs
    DOUBLE PRECISION, DIMENSION(:) :: x

    rows = size(x,1)

    WRITE(*,fmtstr(fs,1)) (x(j), j=1,rows)

  END SUBROUTINE displayvec

! @\newpage\subsection{printmat}@
  SUBROUTINE printmat(x)
    IMPLICIT NONE

    INTENT(IN) :: x
    INTEGER rows,cols,j
    CHARACTER (LEN=10) :: fs
    DOUBLE PRECISION, DIMENSION(:,:) :: x

    fs ='e12.4'
        
    cols = size(x,2)
    rows = size(x,1)

    WRITE(*,fmtstr(fs,cols)) (x(j,:), j=1,rows)

  END SUBROUTINE printmat

! @\newpage\subsection{printvec}@
  SUBROUTINE printvec(x)
    IMPLICIT NONE

    INTENT(IN) :: x
    INTEGER rows,j
    CHARACTER (LEN=200) :: fs
    DOUBLE PRECISION, DIMENSION(:) :: x

    fs ='f20.10'
        
    rows = size(x,1)

    WRITE(*,fmtstr(fs,1)) (x(j), j=1,rows)

  END SUBROUTINE printvec

! @\newpage\subsection{int2str}@
  FUNCTION int2str(n) result(str)
    INTENT(IN) :: n
    INTEGER :: n
    CHARACTER (LEN=200) :: str
    WRITE(str,'(i10)') n
    str = trim(adjustl(str))
  END FUNCTION int2str

! @\newpage\subsection{es30d16}@
  FUNCTION es30d16(N) result(fmtstr)
    INTENT(IN) :: N
    INTEGER :: n
    CHARACTER (LEN=1000) :: fmtstr
    fmtstr = '(ES30.16' // repeat(',ES30.16', N-1) // ')'
  END FUNCTION es30d16

! @\newpage\subsection{nummer84}@
  FUNCTION nummer84(N) result(fmtstr)
    INTENT(IN) :: N
    INTEGER :: n
    CHARACTER (LEN=1000) :: fmtstr
    fmtstr = '(f8.4' // repeat(',f8.4', N-1) // ')'
  END FUNCTION nummer84

! @\newpage\subsection{fmtstr}@
  FUNCTION fmtstr(str,N) 
    INTENT(IN) :: N,str
    INTEGER :: n
    CHARACTER (LEN=1000) :: fmtstr
    CHARACTER (LEN=*) :: str
    fmtstr = '(' // trim(str) // repeat(',' // trim(str), N-1) // ')'
  END FUNCTION fmtstr

! @\newpage\subsection{savemat}@
  SUBROUTINE savemat (x,filename)
    IMPLICIT NONE

    INTEGER :: rows, cols
    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: x(:,:)
    INTEGER j


    rows = size(x,1)
    cols = size(x,2)

    ! WRITE (*,*) 'THIS IS SAVEMAT WITH FILENAME ', filename

    !Open File for writing
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(ES30.16' // repeat(',ES30.16', cols-1) // ')') (x(j,:), j=1,rows)
    CLOSE(UNIT=4)

  END SUBROUTINE savemat

! @\newpage\subsection{savearray3}@
  SUBROUTINE savearray3(x,label,filext)
    IMPLICIT NONE

    INTEGER :: rows, cols, obs
    CHARACTER (LEN=*), INTENT(IN) :: label, filext
    DOUBLE PRECISION, INTENT(IN) :: x(:,:,:)
    CHARACTER (LEN=200) :: filename
    INTEGER jj,ii


    rows = size(x,1)
    cols = size(x,2)
    obs  = size(x,3)

    DO ii=1,obs

       filename = trim(label) // trim(int2str(ii)) // "." // trim(filext)
       
       !Open File for writing
       OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
       WRITE(4,'(ES30.16' // repeat(',ES30.16', cols-1) // ')') (x(jj,:,ii), jj=1,rows)
       CLOSE(UNIT=4)

    END DO
 
  END SUBROUTINE savearray3

! @\newpage\subsection{savematlogical}@
  SUBROUTINE savematlogical(x,filename)
    IMPLICIT NONE

    INTEGER :: rows, cols
    CHARACTER (LEN=*), INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: x(:,:)
    INTEGER j

    rows = size(x,1)
    cols = size(x,2)

    !Open File for writing
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(I5' // repeat(',I5', cols-1) // ')') (merge(1,0,x(j,:)), j=1,rows)
    CLOSE(UNIT=4)

  END SUBROUTINE savematlogical

! @\newpage\subsection{savevec}@
  SUBROUTINE savevec (y,filename)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: y(:)
    INTEGER j

    !Open File for writing
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    DO j=1,size(y)
       WRITE(4,'(ES30.16)') y(j)
    END DO
    CLOSE(UNIT=4)

  END SUBROUTINE savevec

! @\newpage\subsection{savevecX}@
  SUBROUTINE savevecX (y,filename)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: y(:)
    INTEGER j

    !Open File for writing
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    DO j=1,size(y)
       WRITE(4,'(ES40.16E3)') y(j)
    END DO
    CLOSE(UNIT=4)

  END SUBROUTINE savevecX

! @\newpage\subsection{savescalar}@
  SUBROUTINE savescalar (y,filename)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: y

    !Open File for writing
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(ES30.16)') y
    CLOSE(UNIT=4)

  END SUBROUTINE savescalar

! @\newpage\subsection{storeEstimates}@
  SUBROUTINE storeEstimates(theta,Ntheta,Ndraws,filename)
    ! store mean, median and quantiles of draws into file

    IMPLICIT NONE

    INTENT(INOUT) :: theta ! has to be INOUT because of sorting
    INTENT(IN)    :: Ntheta,Ndraws,filename

    INTEGER, parameter :: Nfrac = 10
    REAL, DIMENSION(Nfrac), PARAMETER :: fractiles = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)
    INTEGER, DIMENSION(Nfrac) :: fracndx 
    INTEGER :: Ntheta, Ndraws, n, k, j, status
    DOUBLE PRECISION, DIMENSION(Ntheta,Ndraws) :: theta
    DOUBLE PRECISION, DIMENSION(Ntheta) :: mittel, median

    CHARACTER (LEN=200) :: filename

    ! sort draws
    DO n = 1, Ntheta
       CALL dlasrt('I', Ndraws, theta(n,:), status)
       IF (status /= 0) THEN
          write (*,*), 'DLASORT ERROR ', status, ' [STORE ESTIMATES]'
          stop 1
       END IF
    END DO

    ! generate index for each fractile
    FORALL (n = 1:Nfrac)
       fracndx(n)  = floor(real(Ndraws) * fractiles(n))
    END FORALL
    where (fracndx == 0) fracndx = 1

    ! compute results
    mittel    = sum(theta, 2) / dble(Ndraws)
    median    = theta(:,floor(real(Ndraws) * 0.5))

    ! write to file
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16)') (mittel(j), median(j), (theta(j,fracndx(k)), k = 1,Nfrac), j=1,Ntheta)
    CLOSE(UNIT=4)

  END SUBROUTINE storeEstimates
  ! -----------------------------------------------------------------
  ! @\newpage\subsection{storeEstimatesOMP}@
  SUBROUTINE storeEstimatesOMP(theta,Ntheta,Ndraws,filename)
    ! store mean, median and quantiles of draws into file

    IMPLICIT NONE

    INTENT(INOUT) :: theta ! has to be INOUT because of sorting
    INTENT(IN)    :: Ntheta,Ndraws,filename

    INTEGER, parameter :: Nfrac = 10
    REAL, DIMENSION(Nfrac), PARAMETER :: fractiles = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)
    INTEGER, DIMENSION(Nfrac) :: fracndx 
    INTEGER :: Ntheta, Ndraws, n, k, j, status
    DOUBLE PRECISION, DIMENSION(Ntheta,Ndraws) :: theta
    DOUBLE PRECISION, DIMENSION(Ntheta) :: mittel, median

    CHARACTER (LEN=200) :: filename

    ! sort draws
    !$OMP PARALLEL DO SHARED(Ndraws,theta) PRIVATE(status)
    DO n = 1, Ntheta
       CALL dlasrt('I', Ndraws, theta(n,:), status)
       IF (status /= 0) THEN
          write (*,*), 'DLASORT ERROR ', status, ' [STORE ESTIMATES]'
          stop 1
       END IF
    END DO
    !$OMP END PARALLEL DO

    ! generate index for each fractile
    FORALL (n = 1:Nfrac)
       fracndx(n)  = floor(real(Ndraws) * fractiles(n))
    END FORALL

    ! compute results
    mittel    = sum(theta, 2) / dble(Ndraws)
    median    = theta(:,floor(real(Ndraws) * 0.5))

    ! write to file
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16)') (mittel(j), median(j), (theta(j,fracndx(k)), k = 1,Nfrac), j=1,Ntheta)
    CLOSE(UNIT=4)

  END SUBROUTINE storeEstimatesOMP
  ! -----------------------------------------------------------------
END MODULE embox

