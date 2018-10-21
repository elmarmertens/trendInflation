MODULE blaspack

  use embox , only : savemat, savevec, pi 

  IMPLICIT NONE

CONTAINS


! @\newpage\subsection{symmetric}@
  PURE SUBROUTINE symmetric(S)
    ! ensures symmetry of S (assumung S has at least upper storage)
    IMPLICIT NONE
    INTENT(INOUT) :: S

    INTEGER :: j,i,N
    DOUBLE PRECISION, DIMENSION(:,:) :: S

    N = size(S,1)

       DO i=1,N
          DO j=1,i-1
             S(i,j) = S(j,i)
          END DO
       END DO
  END SUBROUTINE SYMMETRIC

! @\newpage\subsection{maxroot}@
  FUNCTION maxroot(A,n)
    INTENT(IN) :: A, n

    double precision :: maxroot

    integer :: N, status, lwork
    double precision, dimension(n,n) :: A, Awork
    double precision, dimension(n) :: lambdaR, lambdaI

    double precision :: dummy(1,1)
    double precision, allocatable, dimension(:) :: WORK

    Awork = A ! do not touch A 

    ! workspace query
    ALLOCATE (WORK(1))

    LWORK = -1
    call DGEEV('N', 'N', n , Awork, n, lambdaR, lambdaI, dummy, 1, dummy, 1, WORK, LWORK, status)
    IF (status /= 0) THEN
       WRITE(*,*) 'DGEEV error (LWORK QUERY)'
       STOP 1
    END IF

    LWORK = ceiling(WORK(1))

    DEALLOCATE(WORK)

    ! setup optimal workspace
    ALLOCATE(WORK(LWORK))
    ! compute eigenvalues
    call DGEEV('N', 'N', n , Awork, n, lambdaR, lambdaI, dummy, 1, dummy, 1, WORK, LWORK, status)
    IF (status /= 0) THEN
       WRITE(*,*) 'DGEEV error, status =', status
       print *, 'WORK:', WORK(1:5)
       print *, 'LWORK:', LWORK
       call savemat(Awork, 'debug.A.dat')
       STOP 1
    END IF

    maxroot      = maxval(abs(dcmplx(lambdaR, lambdaI)))
    ! wrap up
    DEALLOCATE(WORK)

  END FUNCTION maxroot
! @\newpage\subsection{sandwich}@
  SUBROUTINE sandwich(ASA, A, LDA, S, LDS)
    ! ASA = A * S * A' 
    ! input S is upper-triangular symmetric 
    ! output ASA is dense-symmetric
    INTENT(OUT) :: ASA
    INTENT(IN)  :: A, LDA, S, LDS

    INTEGER :: LDA, LDS

    DOUBLE PRECISION :: A(LDA, LDS), S(LDS,LDS), ASA(LDA,LDA), AS(LDA,LDS)

    AS  = 0.0d0
    ASA = 0.0d0

    call DSYMM('R','U', LDA, LDS,1.0d0,S,LDS,A,LDA,0.0d0,AS,LDA)
    call DGEMM('N','T',LDA,LDA,LDS,1.0d0,AS,LDA,A,LDA,0.0d0,ASA,LDA)

  END SUBROUTINE sandwich

! @\newpage\subsection{sandwichplus}@
  SUBROUTINE sandwichplus(ASA, A, LDA, S, LDS)
    ! ASA = A * S * A' + ASA where S is symmetric (uppper triangular)
    ! input S is upper-triangular symmetric 
    ! input/output ASA has to be dense-symmetric !!

    INTENT(INOUT) :: ASA
    INTENT(IN)  :: A, LDA, S, LDS

    INTEGER :: LDA, LDS

    DOUBLE PRECISION :: A(LDA, LDS), S(LDS,LDS), ASA(LDA,LDA), AS(LDA,LDS)

    AS  = 0.0d0

    call DSYMM('R','U', LDA, LDS,1.0d0,S,LDS,A,LDA,0.0d0,AS,LDA)
    call DGEMM('N','T',LDA,LDA,LDS,1.0d0,AS,LDA,A,LDA,1.0d0,ASA,LDA)

  END SUBROUTINE sandwichplus

! @\newpage\subsection{XprimeX}@
  SUBROUTINE XprimeX(XX, X)

    INTENT(OUT) :: XX
    INTENT(IN)  :: X

    INTEGER :: N, T

    DOUBLE PRECISION, DIMENSION(:,:) :: XX, X

    N = size(X,2)
    T = size(X,1)
    XX = 0.0d0 ! to clean out lower triangular part of XX
    call DSYRK('U','T',N,T,1.0d0,X,T,0.0d0,XX,N)

  END SUBROUTINE XprimeX

! @\newpage\subsection{XXprime}@
  SUBROUTINE XXprime(XX, X)

    INTENT(OUT) :: XX
    INTENT(IN)  :: X

    INTEGER :: Ncols, Nrows

    DOUBLE PRECISION, DIMENSION(:,:) :: XX, X

    Ncols = size(X,2)
    Nrows = size(X,1)
    XX = 0.0d0 ! to clean out lower triangular part of XX
    call DSYRK('U','N',Nrows,Ncols,1.0d0,X,Nrows,0.0d0,XX,Nrows)

  END SUBROUTINE XXprime


! @\newpage\subsection{invsym}@
  SUBROUTINE invsym(xx)
    ! inverts p.d. symmetric real matrix, assuming upper triangular storage

    INTENT(INOUT) :: XX

    INTEGER :: n, info
    DOUBLE PRECISION, DIMENSION(:,:) :: XX

    n = size(xx,1)
    call DPOTRF('U', n, XX, n, info )
    IF (info /= 0) THEN
       write(*,*) "DPOTRF ERROR:", INFO, "[INVSYM]"
       STOP 1
    END IF
    
    call DPOTRI('U', n, XX, n, info )
    IF (info /= 0) THEN
       write(*,*) "DPOTRI ERROR:", INFO, "[INVSYM]"
       STOP 1
    END IF

  END SUBROUTINE invsym


! @\newpage\subsection{symkronecker}@
  SUBROUTINE symkronecker(alpha,A,Na,B,Nb,beta,C)
    ! C = kron(A,B) + C
    ! assumes symmetry A,B and C 
    ! notice: A can be upper triangular, but B must be full storage

    INTENT(OUT) :: C
    INTENT(IN) :: A,Na,B,Nb,alpha,beta

    INTEGER :: Na,Nb, i, j

    DOUBLE PRECISION :: A(Na,Na), B(Nb,Nb), C(Na * Nb, Na * Nb), alpha, beta

    ! loop over rows and columns of A
    DO j = 1 , Na
       FORALL(i=1:j) C((i-1) * Nb + 1 : i * Nb, (j-1) * Nb + 1 : j * Nb) = alpha * A(i,j) * B + beta * C((i-1) * Nb + 1 : i * Nb, (j-1) * Nb + 1 : j * Nb)
    END DO


  END SUBROUTINE symkronecker

! @\newpage\subsection{eye}@
  SUBROUTINE eye(I,alpha)
    ! identity matrix of order n, scaled by alpha (default = 1.0d0)

    INTEGER :: N
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(:,:) :: I
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: alpha
    DOUBLE PRECISION :: a
    INTEGER :: ii


    IF (.NOT. PRESENT(alpha)) THEN
       a = 1.0d0
    ELSE
       a = alpha 
    END IF

    N = size(I,1)
    I = 0.0d0
    FORALL (ii=1:N) I(ii,ii) = a

  END SUBROUTINE eye

! @\newpage\subsection{vec}@
  SUBROUTINE vec(v,x)
    ! v = vec(x)

    INTENT(OUT) :: v
    INTENT(IN) :: x

    INTEGER :: i,j,rows,cols

    DOUBLE PRECISION, DIMENSION(:,:) :: x
    DOUBLE PRECISION, DIMENSION(:) :: v

    rows = size(x,1)
    cols = size(x,2)

    FORALL (i=1:rows,j=1:cols) v((j-1) * rows + i) = x(i,j)

  END SUBROUTINE vec

! @\newpage\subsection{vech}@
  SUBROUTINE vech(v,x)
    ! v = vech(x)
    ! assuming upper triangular storage 

    INTENT(INOUT) :: v
    INTENT(IN) :: x

    INTEGER :: i,j,n,s

    DOUBLE PRECISION, DIMENSION(:,:) :: x
    DOUBLE PRECISION, DIMENSION(:) :: v

    n = size(x,1)
    s = 0
    
    DO j=1,n
       DO i = 1, j
          s = s + 1
          v(s) = x(i,j)
       END DO
    END DO

  END SUBROUTINE vech
! @\newpage\subsection{triu}@
  pure SUBROUTINE triu(x)
    ! zeros out lower triangular elements

    INTENT(INOUT) :: x

    INTEGER :: row,col,n

    DOUBLE PRECISION, DIMENSION(:,:) :: x

    n = size(x,1)

    DO col=1,n-1
       DO row = col+1,n
          x(row,col) = 0.0d0
       END DO
    END DO

  END SUBROUTINE triu

! @\newpage\subsection{choleski}@
  SUBROUTINE choleski(s)
    intent(inout) :: s

    integer :: n, i
    double precision :: s(:,:)

    n = size(s,1)

    ! factorize
    call dpotrf('u', n, s, n, i)

    ! check for errors
    if (i /= 0) then
       write(*,*) 'CHOLESKI ERROR:', i, '[CHOLESKI BLASPACKBOX]'
       stop 1
    end if

    ! zero out lower triangular
    forall (i = 1 : n-1) s(i+1:n,i) = 0.0d0


  END SUBROUTINE choleski

END MODULE blaspack
