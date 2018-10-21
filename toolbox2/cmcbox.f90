MODULE cmcbox

  IMPLICIT NONE

CONTAINS

! @\newpage\subsection{storeEstimatesCMC}@
SUBROUTINE storeEstimatesCMC(theta,Ntheta,Ndraws,filename)
  ! store mean, median and quantiles of draws into file

  IMPLICIT NONE

  INTENT(INOUT) :: theta ! has to be INOUT because of sorting
  INTENT(IN)    :: Ntheta,Ndraws,filename

  INTEGER, parameter :: Nfrac = 9
  REAL, DIMENSION(Nfrac), PARAMETER :: fractiles = (/ 0.025, 0.05, .16, .25, 0.5, .75, .84, 0.95, 0.975 /) ! as obtained from Luke

  INTEGER, DIMENSION(Nfrac) :: fracndx 
  INTEGER :: Ntheta, Ndraws, n, k, j, status
  DOUBLE PRECISION, DIMENSION(Ntheta,Ndraws) :: theta
  DOUBLE PRECISION, DIMENSION(Ntheta) :: mittel !, median

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
  ! median    = theta(:,floor(real(Ndraws) * 0.5))

  ! write to file
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16)') (mittel(j), (theta(j,fracndx(k)), k = 1,Nfrac), j=1,Ntheta)
  CLOSE(UNIT=4)

END SUBROUTINE storeEstimatesCMC
! -----------------------------------------------------------------

END MODULE cmcbox


