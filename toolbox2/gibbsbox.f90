MODULE gibbsbox

  USE blaspack, only : eye, vec, symmetric, symkronecker, XprimeX, xxprime, xxprime, invsym, choleski, maxroot
  USE statespacebox, only : samplerA3B3C3noise, dlyap
  USE embox, only : savemat, savevec, savearray3 ! debugging
  USE timerbox
  USE vslbox

  IMPLICIT NONE

CONTAINS

! @\newpage\subsection{drawRW}@
  SUBROUTINE drawRW (h,T,sigma,Eh0,Vh0,VSLstream)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: T
    DOUBLE PRECISION, DIMENSION(0:T), INTENT(OUT) :: h
    DOUBLE PRECISION, INTENT(IN)  :: sigma, Eh0, Vh0
    DOUBLE PRECISION, DIMENSION(0:T) :: z
    TYPE (vsl_stream_state), INTENT(INOUT) :: VSLstream
    INTEGER :: j, errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T + 1, z, 0.0d0, sigma)

    h(0)    = Eh0 + sqrt(Vh0) / sigma * z(0)
    DO j=1,T
       h(j) = h(j-1) + z(j)
    END DO

  END SUBROUTINE drawRW
! @\newpage\subsection{drawRWcorrelated}@
  SUBROUTINE drawRWcorrelated (h,N,T,sqrtVh,Eh0,sqrtVh0,VSLstream)
    IMPLICIT NONE

    INTENT(IN) :: N,T,sqrtVh,Eh0,sqrtVh0
    INTENT(OUT) :: h
    INTENT(INOUT) :: VSLstream

    INTEGER :: T,N, j, errcode
    DOUBLE PRECISION, DIMENSION(N,0:T) :: h, z
    DOUBLE PRECISION :: sqrtVh(N,N), Eh0(N), sqrtVh0(N,N)
    TYPE (vsl_stream_state) :: VSLstream
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

    ! draw random numbers
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * N, z, 0.0d0, 1.0d0)
    ! construct h0
    h(:,0) = Eh0
    call DGEMV('N',N,N,1.0d0,sqrtVh0,N,z(:,0),1,1.0d0,h(:,0),1)
    ! scale shocks and accumulate
    DO j=1,T
       h(:,j) = h(:,j-1) 
       call DGEMV('N',N,N,1.0d0,sqrtVh,N,z(:,j),1,1.0d0,h(:,j),1)
    END DO

  END SUBROUTINE drawRWcorrelated
! @\newpage\subsection{smoothingsamplerLocalLevel}@
  SUBROUTINE smoothingsamplerLocalLevel (tau,y,T,vartrend,varnoise,X0,V0, VSLstream)
    IMPLICIT NONE
    INTENT (IN) :: T, y, vartrend, varnoise, X0, V0
    INTENT (OUT) :: tau
    INTENT(INOUT) :: VSLstream
    INTEGER j, errcode, T
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    DOUBLE PRECISION, DIMENSION(0:T) :: tautT, tau, tauplus, SigmaStar, z
    DOUBLE PRECISION, DIMENSION(T) :: y, yplus, vartrend, varnoise, Sigma, e
    DOUBLE PRECISION :: X0, V0, gain
    TYPE (vsl_stream_state)  :: VSLstream


    ! draw plus data
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T+1, z, 0.0d0, 1.0d0)
    tauplus(0) = X0 + sqrt(V0) * z(0)
    DO j=1,T
       tauplus(j) = tauplus(j-1) + sqrt(vartrend(j)) * z(j)
    END DO
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T, e, 0.0d0, 1.0d0)

    ! Prepare filter for "pseudo" observables
    yplus = tauplus(1:T) + sqrt(varnoise) * e
    yplus = yplus - y

    SigmaStar(0)  = V0
    tautT(0)      = 0.0d0 ! adjust for yplus minus y

    ! forward filter
    DO j=1,T

       Sigma(j)       = SigmaStar(j-1) + vartrend(j)
       gain           = Sigma(j) / (Sigma(j) + varnoise(j))

       SigmaStar(j)   = ((1 - gain)**2) * Sigma(j) + (gain**2) * varnoise(j)

       tautT(j)       = (1 - gain) * tautT(j-1) + gain * yplus(j)

    END DO

    ! backward filter
    DO j=T-1,0,-1 
       gain     = SigmaStar(j) / Sigma(j+1)
       tautT(j) = (1 - gain) * tautT(j) + gain * tautT(j+1)
    END DO

    ! put draws together
    tau = tauplus - tautT

    ! DEBUG
    ! OPEN (UNIT=4, FILE='debug.smoothingsamplerSettings.dat', STATUS='REPLACE', ACTION='WRITE')
    ! WRITE(4,'(2ES30.16,I10)') X0,V0,T
    ! CLOSE(UNIT=4)
    ! OPEN (UNIT=4, FILE='debug.smoothingsamplerT.dat', STATUS='REPLACE', ACTION='WRITE')
    ! WRITE(4,'(6ES30.16)') (y(j), yplus(j), varnoise(j), vartrend(j), Sigma(j), e(j), j=1,T)
    ! CLOSE(UNIT=4)
    ! OPEN (UNIT=4, FILE='debug.smoothingsamplerTp1.dat', STATUS='REPLACE', ACTION='WRITE')
    ! WRITE(4,'(5ES30.16)') (tau(j), tauplus(j), tautT(j), SigmaStar(j), z(j), j=0,T)
    ! CLOSE(UNIT=4)

  END SUBROUTINE smoothingsamplerLocalLevel

! @\newpage\subsection{tvpRegressionSlope}@
  SUBROUTINE tvpRegressionSlope(beta,T,y,x,varbeta,varnoise,beta0,beta0V,VSLstream)
    ! univarate regression with tvp slope 
    ! y(t) = beta(t) * x(t) + sqrt(varnoise(t)) * e(t)
    ! beta(t) = beta(t-1) + sqrt(varbeta(t)) * z(t)

    IMPLICIT NONE
    INTENT (IN) :: T, y, x, varbeta, varnoise, beta0, beta0V
    INTENT (OUT) :: beta
    INTENT(INOUT) :: VSLstream
    INTEGER j, errcode, T
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    DOUBLE PRECISION, DIMENSION(0:T) :: betatT, beta, betaplus, SigmaStar, z
    DOUBLE PRECISION, DIMENSION(T) :: y, yplus, x, Sigma, e
    DOUBLE PRECISION :: beta0, beta0V, gain, omxgain, varbeta, varnoise
    TYPE (vsl_stream_state)  :: VSLstream


    ! draw plus data
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T+1, z, 0.0d0, sqrt(varbeta))
    betaplus(0) = beta0 + sqrt(beta0V  / varbeta) * z(0)
    DO j=1,T
       betaplus(j) = betaplus(j-1) + z(j)
    END DO

    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T, e, 0.0d0, sqrt(varnoise))
    ! Prepare filter for "pseudo" observables
    yplus = x * betaplus(1:T) + e
    yplus = yplus - y

    SigmaStar(0)  = beta0V
    betatT(0)     = 0.0d0 ! adjust for yplus minus y

    ! forward filter
    DO j=1,T

       Sigma(j)       = SigmaStar(j-1) + varbeta
       gain           = x(j) * Sigma(j) / (Sigma(j) * (x(j) ** 2) + varnoise)
       omxgain        = 1 - gain * x(j)

       SigmaStar(j)   = omxgain * Sigma(j) 
       betatT(j)      = omxgain * betatT(j-1) + gain * yplus(j)

    END DO

    ! backward filter
    DO j=T-1,0,-1 
       gain     = SigmaStar(j) / Sigma(j+1)
       betatT(j) = (1 - gain) * betatT(j) + gain * betatT(j+1)
    END DO

    ! put draws together
    beta = betaplus - betatT

    ! DEBUG
    ! call savevec(beta, 'tvpdebug.beta.dat')
    ! call savevec(y, 'tvpdebug.y.dat')
    ! call savevec(x, 'tvpdebug.x.dat')
    ! call savevec(varbeta, 'tvpdebug.varbeta.dat')
    ! call savevec(varnoise, 'tvpdebug.varnoise.dat')
    ! call savevec(z, 'tvpdebug.z.dat')
    ! call savevec(e, 'tvpdebug.e.dat')

  END SUBROUTINE tvpRegressionSlope

! @\newpage\subsection{GelmanTest1}@
  SUBROUTINE GelmanTest1(SRstat,draws,Nsims,Nstreams)
    IMPLICIT NONE

    INTENT(OUT) :: SRstat
    INTENT(IN) :: Nsims, Nstreams
    INTENT(IN) :: draws

    INTEGER :: Nsims, Nstreams
    DOUBLE PRECISION, DIMENSION(Nsims,Nstreams)  :: draws
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: edraws
    DOUBLE PRECISION :: SRstat, W, B, psibar
    DOUBLE PRECISION, DIMENSION(Nstreams) :: psi, s2 ! mean and variance within stream
    INTEGER :: j,status

    ! within stream means
    psi    = sum(draws,1) / Nsims
    psibar = sum(psi) / Nstreams

    ! variance between means
    B  = sum((psi - psibar) ** 2) / (Nstreams - 1)


    ALLOCATE (edraws(Nsims,Nstreams), STAT=status)
    IF (status /=0) THEN
       WRITE (*,*) 'Could not allocate edraws'
    END IF

    ! within stream variances
    FORALL (j=1:Nstreams)
       edraws(:,j) = draws(:,j) - psi(j)
    END FORALL
    s2 = sum(edraws ** 2, 1) / (Nsims - 1)
    ! average of within-stream-variances
    W  = sum(s2) / Nstreams

    DEALLOCATE (edraws, STAT=status)
    IF (status /=0) THEN
       WRITE (*,*) 'Could not de-allocate edraws'
    END IF

    SRstat = sqrt(dble(Nsims - 1) / dble(Nsims)  + B / W)

  END SUBROUTINE GelmanTest1


! @\newpage\subsection{drawNDXpdf0}@
  SUBROUTINE drawNDXpdf0(ndx, Ndraws, pdf, N, VSLstream)
    ! draw indices of multinominal distribution (given pdf)
    ! NOTE: the implementation via COUNT seems to be a little faster, since it needs to compute the cdf only once
    IMPLICIT NONE
    INTENT(OUT) :: ndx
    INTENT(IN) :: Ndraws, N, pdf
    INTENT(INOUT) :: VSLstream

    INTEGER :: Ndraws, N, errcode, j, i
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    DOUBLE PRECISION, DIMENSION(Ndraws) :: u
    INTEGER, DIMENSION(Ndraws) :: ndx
    DOUBLE PRECISION :: pdf(N), cdf ! , cdf2(N)
    TYPE (vsl_stream_state) :: VSLstream

    ! draw uniforms
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, Ndraws, u, 0.0d0, 1.0d0)

    ndx = 0
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(u, ndx, pdf, Ndraws) PRIVATE(cdf,i)
    DO j=1,Ndraws
       cdf = 0.0d0
       i   = 0
       DO WHILE (cdf < u(j))
          i      = i + 1
          cdf    = cdf + pdf(i)
       END DO
       ndx(j) = i
    END DO
    !$OMP END PARALLEL DO 

  END SUBROUTINE drawNDXpdf0

! @\newpage\subsection{drawNDXpdf}@
  SUBROUTINE drawNDXpdf(ndx, Ndraws, pdf, N, VSLstream)
    ! draw indices of multinominal distribution (given pdf)
    IMPLICIT NONE
    INTENT(OUT) :: ndx
    INTENT(IN) :: Ndraws, N, pdf
    INTENT(INOUT) :: VSLstream

    INTEGER :: Ndraws, N, errcode, j
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    DOUBLE PRECISION, DIMENSION(Ndraws) :: u
    INTEGER, DIMENSION(Ndraws) :: ndx
    DOUBLE PRECISION, DIMENSION(N) :: pdf, cdf
    TYPE (vsl_stream_state) :: VSLstream


    ! cumulate pdf into cdf
    cdf(1) = pdf(1)
    DO j=2,N
       cdf(j) = pdf(j) + cdf(j-1)
    END DO

    ! draw uniforms
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, Ndraws, u, 0.0d0, 1.0d0)

    ! sort the uniforms (in order to obtain sorted ndx)
    ! CALL dlasrt('I', Ndraws, u, errcode)
   
    !$OMP PARALLEL DO SHARED(ndx, u, cdf)
    DO j=1,Ndraws
       ndx(j) = COUNT(u(j) > cdf) + 1
    END DO
    !$OMP END PARALLEL DO 

  END SUBROUTINE drawNDXpdf

! @\newpage\subsection{drawNDXpdf2}@
! @\newpage\subsection{drawKSCstates}@
  SUBROUTINE drawKSCstates(kai2states, T, N, VSLstream)
    ! unconditional states draws of KSC7 mixture probability
    ! this is a vectorized version of drawKai2statesKSC7
    IMPLICIT NONE
    INTENT(OUT) :: kai2states
    INTENT(IN) :: T, N
    INTENT(INOUT) :: VSLstream

    INTEGER :: T, N, errcode, j, i
    INTEGER, PARAMETER :: KSCmix = 7, VSLmethodGaussian = 0, VSLmethodUniform = 0
    DOUBLE PRECISION, DIMENSION(T,N) :: u
    INTEGER, DIMENSION(T,N) :: kai2states
    DOUBLE PRECISION, DIMENSION(KSCmix) :: KSCcdf
    TYPE (vsl_stream_state) :: VSLstream

 
    KSCcdf = (/ 7.3d-3, 112.86d-3, 112.88d-3, 156.83d-3, 496.84d-3, 742.5d-3, 1.0d0 /)

    ! KSCcdf     = (/ .0073d0, .10556d0, .00002d0, .04395d0, .34001d0, .24566d0, .25750d0 /)
    ! DO j=2,KSCmix
    !    KSCcdf(j) = KSCcdf(j) + KSCcdf(j-1)
    ! END DO

    errcode = vdrnguniform(VSLmethodUniform, VSLstream, T * N, u, 0.0d0, 1.0d0)
    FORALL (j=1:T, i =1:N)
       kai2states(j,i) = COUNT(u(j,i) > KSCcdf) + 1
    END FORALL

  END SUBROUTINE drawKSCstates
! @\newpage\subsection{stochvolKSC0}@
  SUBROUTINE stochvolKSC0(h, y, T, hInno, Eh0, Vh0, VSLstream)
    
    ! uses corrected MCMC order as per DelNegro and Primiceri (2013)
    ! same as stochvolKSC, except that kai2states are not used/stored as argument anymore
    ! ... and slight change in the order of arguments ...

    IMPLICIT NONE

    INTENT(IN) :: y, hInno, Eh0, Vh0, T
    INTENT(INOUT) :: VSLstream
    INTENT(INOUT) :: h

    INTEGER :: T, errcode, j, s
    INTEGER, PARAMETER :: KSCmix = 7
    DOUBLE PRECISION, DIMENSION(KSCmix), PARAMETER :: KSCmean = - 1.2704d0 + (/ -10.12999d0, -3.97281d0, -8.56686d0, 2.77786d0, .61942d0, 1.79518d0, -1.08819d0 /), KSCvar     = (/ 5.79596d0, 2.61369d0, 5.1795d0, .16735d0, .64009d0, .34023d0, 1.26261d0 /), KSCpdf= (/ .0073d0, .10556d0, .00002d0, .04395d0, .34001d0, .24566d0, .25750d0 /), KSCvol = sqrt(KSCvar)

    DOUBLE PRECISION, DIMENSION(T) :: y, logy2, logy2star, varnoisemix, varh, u
    DOUBLE PRECISION, DIMENSION(T,KSCmix) :: kai2CDF
    INTEGER, DIMENSION(T) :: kai2states
    DOUBLE PRECISION, DIMENSION(0:T) :: h
    DOUBLE PRECISION  :: hInno, Eh0, Vh0

    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream

    ! log-linear observer
    logy2 = log(y ** 2 + 0.001d0)

    ! STEP 1 DRAW KAI2STATES
    ! a) construct PDF for draws (stored in kai2CDF)
    FORALL (s=1:KSCmix,j=1:T)
       kai2CDF(j,s) = exp(-0.5d0 * ((logy2(j) - h(j) - KSCmean(s)) / KSCvol(s))** 2) / KSCvol(s) * KSCpdf(s)
    END FORALL

    ! b) convert PDF into CDF for draws
    DO s=2,KSCmix
       kai2CDF(:,s) = kai2CDF(:,s-1) + kai2CDF(:,s)
    END DO
    DO s=1,KSCmix-1
       kai2CDF(:,s) = kai2CDF(:,s) / kai2CDF(:,KSCmix)
    END DO
    kai2CDF(:,KSCmix) = 1.0d0


    ! c) draw kai2states
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, T, u, 0.0d0, 1.0d0 )
    FORALL (j=1:T)
       kai2states(j) = COUNT(u(j) > kai2CDF(j,:)) + 1
    END FORALL


    ! STEP 2: KALMAN FILTER FOR h
    
    ! construct trend variance and noise variance
    varh = hInno ** 2
    FORALL (j=1:T)
       varnoisemix(j) = KSCvar(kai2states(j))
    END FORALL

    ! demeaned observables 
    FORALL(j=1:T)
       logy2star(j) = logy2(j) - KSCmean(kai2states(j))
    END FORALL
    CALL smoothingsamplerLocalLevel(h,logy2star,T,varh,varnoisemix,Eh0,Vh0,VSLstream)

  END SUBROUTINE stochvolKSC0

! @\newpage\subsection{stochvolKSCjprwrap}@
  SUBROUTINE stochvolKSCjprwrap(SVol, h, y, T, hInno, Eh0, Vh0, VSLstream)
    
    ! uses corrected MCMC order as per DelNegro and Primiceri (2013)
    ! same as stochvolKSC, except that kai2states are not used/stored as argument anymore
    ! ... and slight change in the order of arguments ...
    ! this version works with same inputs/outputs as jpr0

    IMPLICIT NONE

    INTENT(IN) :: y, hInno, Eh0, Vh0, T
    INTENT(INOUT) :: VSLstream, SVol
    INTENT(OUT) :: h

    INTEGER :: T, errcode, j, s
    INTEGER, PARAMETER :: KSCmix = 7
    DOUBLE PRECISION, DIMENSION(KSCmix), PARAMETER :: KSCmean = - 1.2704d0 + (/ -10.12999d0, -3.97281d0, -8.56686d0, 2.77786d0, .61942d0, 1.79518d0, -1.08819d0 /), KSCvar     = (/ 5.79596d0, 2.61369d0, 5.1795d0, .16735d0, .64009d0, .34023d0, 1.26261d0 /), KSCpdf= (/ .0073d0, .10556d0, .00002d0, .04395d0, .34001d0, .24566d0, .25750d0 /), KSCvol = sqrt(KSCvar)

    DOUBLE PRECISION, DIMENSION(T) :: y, logy2, logy2star, varnoisemix, varh, u
    DOUBLE PRECISION, DIMENSION(T,KSCmix) :: kai2CDF
    INTEGER, DIMENSION(T) :: kai2states
    DOUBLE PRECISION, DIMENSION(0:T) :: h, SVol
    DOUBLE PRECISION  :: hInno, Eh0, Vh0

    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream

    ! log-linear observer
    h = 2.0d0 * log(SVol)

    logy2 = log(y ** 2 + 0.001d0)

    ! STEP 1 DRAW KAI2STATES
    ! a) construct PDF for draws (stored in kai2CDF)
    FORALL (s=1:KSCmix,j=1:T)
       kai2CDF(j,s) = exp(-0.5d0 * ((logy2(j) - h(j) - KSCmean(s)) / KSCvol(s))** 2) / KSCvol(s) * KSCpdf(s)
    END FORALL

    ! b) convert PDF into CDF for draws
    DO s=2,KSCmix
       kai2CDF(:,s) = kai2CDF(:,s-1) + kai2CDF(:,s)
    END DO
    DO s=1,KSCmix-1
       kai2CDF(:,s) = kai2CDF(:,s) / kai2CDF(:,KSCmix)
    END DO
    kai2CDF(:,KSCmix) = 1.0d0


    ! c) draw kai2states
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, T, u, 0.0d0, 1.0d0 )
    FORALL (j=1:T)
       kai2states(j) = COUNT(u(j) > kai2CDF(j,:)) + 1
    END FORALL


    ! STEP 2: KALMAN FILTER FOR h
    
    ! construct trend variance and noise variance
    varh = hInno ** 2
    FORALL (j=1:T)
       varnoisemix(j) = KSCvar(kai2states(j))
    END FORALL

    ! demeaned observables 
    FORALL(j=1:T)
       logy2star(j) = logy2(j) - KSCmean(kai2states(j))
    END FORALL
    CALL smoothingsamplerLocalLevel(h,logy2star,T,varh,varnoisemix,Eh0,Vh0,VSLstream)

    SVol = exp(h * 0.5d0)

  END SUBROUTINE stochvolKSCjprwrap

! @\newpage\subsection{igammaDraw}@
  SUBROUTINE igammaDraw(draw, sigma0T, dof0, VSLstream)

    INTENT(INOUT) :: draw, VSLstream
    INTENT(IN)    :: sigma0T, dof0

    INTEGER :: dof0, errcode
    type (vsl_stream_state) :: VSLstream

    DOUBLE PRECISION, DIMENSION(dof0) :: z
    DOUBLE PRECISION :: sigma0T, draw

    INTEGER, PARAMETER :: VSLmethodGaussian = 0

    ! draw
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, dof0, z, 0.0d0, 1.0d0)
    draw = sigma0T / sum(z ** 2)

  END SUBROUTINE igammaDraw

! @\newpage\subsection{varianceDraw}@
  SUBROUTINE varianceDraw(igammadraw, sigma0T, dof0, resid, Nobs, VSLstream)

    INTENT(INOUT) :: igammadraw, VSLstream
    INTENT(IN)    :: resid, Nobs, sigma0T, dof0

    INTEGER :: Nobs, dof, dof0, errcode
    type (vsl_stream_state) :: VSLstream

    DOUBLE PRECISION, DIMENSION(Nobs) :: resid
    DOUBLE PRECISION, DIMENSION(dof0 + Nobs) :: z
    DOUBLE PRECISION :: sigmaT, sigma0T, igammadraw

    INTEGER, PARAMETER :: VSLmethodGaussian = 0

    ! compute posterior
    dof     = dof0 + Nobs
    sigmaT  = sigma0T + sum(resid ** 2)
   
    ! draw
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, dof, z, 0.0d0, 1.0d0)
    igammadraw = sigmaT / sum(z ** 2)

    ! if (igammadraw < 0.0d0) then 
    !    print *, 'houston'
    !    print *, sigma0T
    !    print *, sigmaT
    !    print *, dof
    !    print *, sum(z ** 2)
    !    stop 33
    ! end if
    
  END SUBROUTINE varianceDraw
! @\newpage\subsection{iwishDraw}@
  SUBROUTINE iwishDraw(iwdraw, SigmaT, dof, Ny, VSLstream)

    ! SigmaT is supposed to be in upper triangular storage 

    INTENT(INOUT) :: iwdraw, VSLstream
    INTENT(IN)    :: SigmaT, dof, Ny

    INTEGER :: dof, errcode, Ny
    type (vsl_stream_state) :: VSLstream

    DOUBLE PRECISION, DIMENSION(Ny, dof) :: z
    DOUBLE PRECISION, DIMENSION(Ny,Ny) :: SigmaT, work, iwdraw
    ! INTEGER, DIMENSION(Ny) :: ipiv ! helper variable for dgetrs 

    INTEGER, PARAMETER :: VSLmethodGaussian = 0

    ! ztilde = SigmaT^{-.5}' z
    work = SigmaT
    ! choleski of SigmaT
    call DPOTRF('U', Ny, work, Ny, errcode)
    ! zero out lower triangular portion of the choleski --should not be necessary though
    ! FORALL (ii = 1 : Ny-1) work(ii+1:Ny,ii) = 0 

    ! draw z
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, Ny * dof, z, 0.0d0, 1.0d0)

    ! compute ztilde
    ! FORALL (ii = 1 : Ny) ipiv(ii) = ii
    ! dummy = z
    ! call DGETRS('N', Ny, dof, work, Ny, IPIV, dummy, Ny, errcode) 
    ! call savemat(dummy, 'dummy.debug')
    call DTRTRS('U', 'N', 'N', Ny, dof, work, Ny, z, Ny, errcode)
    IF (errcode /= 0)  THEN
       WRITE(*,*) "DTRTRS error:", errcode, '[iwishdraw]'
       STOP 1
    END IF


    ! wishdraw = ztilde * ztilde'
    call XXprime(iwdraw, z) ! recall: upper triangular storage
    call invsym(iwdraw)

  END SUBROUTINE iwishDraw
! @\newpage\subsection{iwishcholDraw}@
  SUBROUTINE iwishcholDraw(iwcdraw, SigmaT, dof, Ny, VSLstream)
    
    ! returns choleski of iwishdraw (however: chol is LHS-upper triangular)

    ! SigmaT is supposed to be in upper triangular storage 

    INTENT(INOUT) :: iwcdraw, VSLstream
    INTENT(IN)    :: SigmaT, dof, Ny

    INTEGER :: dof, errcode, Ny
    type (vsl_stream_state) :: VSLstream

    DOUBLE PRECISION, DIMENSION(Ny, dof) :: z
    DOUBLE PRECISION, DIMENSION(Ny,Ny) :: SigmaT, iwcdraw
    ! INTEGER, DIMENSION(Ny) :: ipiv ! helper variable for dgetrs 

    INTEGER, PARAMETER :: VSLmethodGaussian = 0

    IF (dof < Ny) THEN
       WRITE(*,*) "IWISHCHOLDRAW: need dof at least as large as Ny"
       STOP 1
    END IF


    ! ztilde = SigmaT^{-.5}' z
    iwcdraw = SigmaT
    ! choleski of SigmaT
    call DPOTRF('U', Ny, iwcdraw, Ny, errcode)
    ! zero out lower triangular portion of the choleski --should not be necessary though
    ! FORALL (ii = 1 : Ny-1) iwcdraw(ii+1:Ny,ii) = 0 

    ! draw z
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, Ny * dof, z, 0.0d0, 1.0d0)

    ! compute ztilde
    ! FORALL (ii = 1 : Ny) ipiv(ii) = ii
    ! call DGETRS('N', Ny, dof, work, Ny, IPIV, z, Ny, errcode) 
    call DTRTRS('U', 'N', 'N', Ny, dof, iwcdraw, Ny, z, Ny, errcode)
    IF (errcode /= 0)  THEN
       WRITE(*,*) "DTRTRS error:", errcode, '[iwishcholdraw]'
       STOP 1
    END IF

    ! work = ztilde * ztilde' (wishart draw)
    iwcdraw = 0.0d0 ! ensures that matrix remains lower triangular
    call DSYRK('U','N',Ny,dof,1.0d0,z,Ny,0.0d0,iwcdraw,Ny)

    ! chol(wishdraw)
    call DPOTRF('U', Ny, iwcdraw, Ny, errcode)
    IF (errcode /= 0) THEN
       write(*,*) "DPOTRF ERROR:", errcode, "[iwishcholDRAW]"
       STOP 1
    END IF

    call DTRTRI('U', 'N', Ny, iwcdraw, Ny, errcode)
    IF (errcode /= 0) THEN
       write(*,*) "DTRTRI ERROR:", errcode, "[iwishcholDRAW]"
       STOP 1
    END IF

    ! call savemat(iwcdraw, 'iwc1.debug')
    ! FORALL (ii = 1 : Ny-1) work(ii+1:Ny,ii) = 0
    ! FORALL (ii = 1 : Ny) ipiv(ii) = ii
    ! call eye(iwcdraw)
    ! call DGETRS ('N', Ny, Ny, work, Ny, IPIV, iwcdraw, Ny,errcode) 
    ! call savemat(iwcdraw, 'iwc2.debug')
    ! stop 22

  END SUBROUTINE iwishcholDraw
! @\newpage\subsection{vcvcholDrawTR}@
  SUBROUTINE vcvcholDrawTR(iwishcholdraw, Sigma0T, dof0, resid, T, Ny, VSLstream)

    ! SigmaT is supposed to be in upper triangular storage

    INTENT(OUT)   :: iwishcholdraw
    INTENT(INOUT) :: VSLstream
    INTENT(IN)    :: resid, T, Ny, Sigma0T, dof0

    INTEGER :: T, dof, dof0, errcode, Ny
    type (vsl_stream_state) :: VSLstream

    DOUBLE PRECISION, DIMENSION(T,Ny) :: resid
    DOUBLE PRECISION, DIMENSION(Ny, dof0 + T) :: z
    DOUBLE PRECISION, DIMENSION(Ny,Ny) :: Sigma0T, SigmaT, iwishcholdraw
    ! INTEGER, DIMENSION(Ny) :: ipiv

    INTEGER, PARAMETER :: VSLmethodGaussian = 0

    ! compute posterior
    dof     = dof0 + T

    ! important: getrs needs upper triangular factorization

    ! SigmaT = resid' * resid + Sigma0T 
    SigmaT = Sigma0T
    call DSYRK('U','T',Ny,T,1.0d0,resid,T,1.0d0,SigmaT,Ny) 

    ! choleski of SigmaT
    call DPOTRF('U', Ny, SigmaT, Ny, errcode)
    IF (errcode /= 0) THEN
       write(*,*) "DPOTRF ERROR:", errcode, "[VCVCHOLDRAW: SIGMAT]"
       STOP 1
    END IF
    ! FORALL (ii = 1 : Ny-1) SigmaT(ii+1:Ny,ii) = 0

    ! ztilde = SigmaT^{-.5}' z
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, Ny * dof, z, 0.0d0, 1.0d0)
    ! FORALL (ii = 1 : Ny) ipiv(ii) = ii
    ! call DGETRS('N', Ny, dof, SigmaT, Ny, IPIV, z, Ny, errcode)
    call DTRTRS('U', 'N', 'N', Ny, dof, SigmaT, Ny, z, Ny, errcode)
    IF (errcode /= 0)  THEN
       WRITE(*,*) "DTRTRS error:", errcode
       STOP 1
    END IF

    ! wishdraw = ztilde * ztilde'
    iwishcholdraw = 0.0d0 ! this ensure that iwishcholdraw is properly lower triangular
    call DSYRK('U','N',Ny,dof,1.0d0,z,Ny,0.0d0,iwishcholdraw,Ny)
    ! chol(wishdraw)
    call DPOTRF('U', Ny, iwishcholdraw, Ny, errcode)
    IF (errcode /= 0) THEN
       write(*,*) "DPOTRF ERROR:", errcode, "[VCVCHOLDRAW: iwishcholdraw]"
       call savemat(resid, 'resid.debug')
       STOP 1
    END IF
    
    ! invert chol(wishdraw)
    call DTRTRI('U', 'N', Ny, iwishcholdraw, Ny, errcode)
    IF (errcode /= 0) THEN
       write(*,*) "DTRTRI ERROR:", errcode, "[VCVCHOLDRAW]"
       STOP 1
    END IF

  END SUBROUTINE vcvcholDrawTR


  ! @\newpage\subsection{SVHdiffusecholeskiKSC}@
  SUBROUTINE SVHdiffusecholeskiKSC(T, Ny, SVol, h, shockslopes, Nslopes, y, hInno, E0h, V0h, VSLstream)

    INTENT(INOUT) :: SVol, VSLstream, y
    INTENT(IN) :: T,Ny, hInno, E0h, V0h
    INTENT(OUT) :: shockslopes, h

    INTEGER :: i, k, T, Ny, Nslopes, offset, these
    DOUBLE PRECISION :: h(Ny,0:T), SVol(Ny,0:T), y(Ny,T), hInno(Ny), E0h(Ny), V0h(Ny), shockslopes(Nslopes)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream


    DOUBLE PRECISION :: lhs(T), rhs(T,Ny-1)

    ! NOTE: V0slopes (and thus also iV0slopes) is assumed to be block-diagonal, with separate blocks for each regression

    ! PART 1: Shockslopes
    offset = 0
    DO i=2,Ny
       these = i - 1

       FORALL (k=1:T) lhs(k)       = y(i,k) / SVol(i,k)
       FORALL (k=1:T) rhs(k,1:i-1) = y(1:i-1,k) / SVol(i,k)

       ! slope indices are offset+1:offset+these
       ! call bayesRegressionSlope(shockslopes(offset+1:offset+these), lhs, rhs(:,1:i-1), these, T, 1.0d0, E0slopes(offset+1:offset+these), iV0slopes(offset+1:offset+these,offset+1:offset+these), VSLstream)
       call bayesdiffuseRegressionSlope(shockslopes(offset+1:offset+these), lhs, rhs(:,1:i-1), these, T, 1.0d0, VSLstream)

       FORALL (k=1:T) y(i,k) = lhs(k) * SVol(i,k)

       offset = offset + these

    END DO
    

    ! PART 2: SV
    DO i = 1,Ny
       CALL stochvolKSCjprwrap(SVol(i,:), h(i,:), y(i,:), T, hInno(i), E0h(i), V0h(i), VSLstream)
    END DO

  END SUBROUTINE SVHdiffusecholeskiKSC

! @\newpage\subsection{SVHcholeskiKSC}@
  SUBROUTINE SVHcholeskiKSC(T, Ny, SVol, h, shockslopes, Nslopes, y, E0slopes, iV0slopes, hInno, E0h, V0h, VSLstream)
    
    ! h = log(SVol ** 2)

    INTENT(INOUT) :: SVol, VSLstream, h
    INTENT(IN) :: T,Ny, E0slopes, iV0slopes, hInno, E0h, V0h, y
    INTENT(OUT) :: shockslopes

    INTEGER :: i, k, T, Ny, Nslopes, offset, these
    DOUBLE PRECISION :: h(Ny,0:T), SVol(Ny,0:T), resid(Ny,T), y(Ny,T), E0slopes(Nslopes), iV0slopes(Nslopes, Nslopes), hInno(Ny), E0h(Ny), V0h(Ny), shockslopes(Nslopes)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream


    DOUBLE PRECISION :: lhs(T), rhs(T,Ny-1)

    ! NOTE: V0slopes (and thus also iV0slopes) is assumed to be block-diagonal, with separate blocks for each regression

    resid = y

    ! PART 1: Shockslopes
    offset = 0
    DO i=2,Ny
       these = i - 1

       FORALL (k=1:T) lhs(k)       = resid(i,k) / SVol(i,k)
       FORALL (k=1:T) rhs(k,1:i-1) = resid(1:i-1,k) / SVol(i,k)

       ! slope indices are offset+1:offset+these
       call bayesRegressionSlope(shockslopes(offset+1:offset+these), lhs, rhs(:,1:i-1), these, T, 1.0d0, E0slopes(offset+1:offset+these), iV0slopes(offset+1:offset+these,offset+1:offset+these), VSLstream)

       FORALL (k=1:T) resid(i,k) = lhs(k) * SVol(i,k)

       offset = offset + these

    END DO
    

    ! PART 2: SV
    DO i = 1,Ny
       CALL stochvolKSCjprwrap(SVol(i,:), h(i,:), resid(i,:), T, hInno(i), E0h(i), V0h(i), VSLstream)
    END DO

  END SUBROUTINE SVHcholeskiKSC

  ! @\newpage\subsection{SVHcholeskiKSCAR1corplus}@
  SUBROUTINE SVHcholeskiKSCAR1corplus(T, Nsv, SVol, h, hshock, hbar, rho, shockslopes, Nslopes, y, E0slopes, iV0slopes, sqrtVhshock, Eh0, sqrtVh0, VSLstream)

    ! h = log(SVol ** 2)

    INTENT(INOUT) :: SVol, VSLstream, y, h
    INTENT(IN) :: T,Nsv, E0slopes, iV0slopes, sqrtVhshock, Eh0, sqrtVh0, rho
    INTENT(OUT) :: shockslopes, hshock, hbar

    DOUBLE PRECISION :: rho 

    INTEGER :: i, s, j, k, T, Nsv, Nslopes, offset, these, errcode
    DOUBLE PRECISION :: E0slopes(Nslopes), iV0slopes(Nslopes, Nslopes), sqrtVhshock(Nsv,Nsv), Eh0(Nsv), sqrtVh0(Nsv,Nsv), shockslopes(Nslopes)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream

    ! KSC mixture
    INTEGER, PARAMETER :: KSCmix = 7
    DOUBLE PRECISION, DIMENSION(KSCmix), PARAMETER :: KSCmean = - 1.2704d0 + (/ -10.12999d0, -3.97281d0, -8.56686d0, 2.77786d0, .61942d0, 1.79518d0, -1.08819d0 /), KSCvar     = (/ 5.79596d0, 2.61369d0, 5.1795d0, .16735d0, .64009d0, .34023d0, 1.26261d0 /), KSCpdf= (/ .0073d0, .10556d0, .00002d0, .04395d0, .34001d0, .24566d0, .25750d0 /), KSCvol = sqrt(KSCvar)
    DOUBLE PRECISION, DIMENSION(Nsv,T,KSCmix) :: kai2CDF
    INTEGER, DIMENSION(Nsv,T) :: kai2states


    DOUBLE PRECISION, DIMENSION(Nsv,T) :: y, logy2, logy2star, volnoisemix, u
    DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: h, SVol
    DOUBLE PRECISION, DIMENSION(Nsv,1:T) :: hshock
    DOUBLE PRECISION, DIMENSION(Nsv) :: hbar

    ! state space matrices
    DOUBLE PRECISION :: A(Nsv*2,Nsv*2,T), B(Nsv*2,Nsv,T), C(Nsv,Nsv*2,T), y2noise(Nsv,T), State(Nsv*2,0:T), StateShock(Nsv*2,1:T), sqrtState0V(2*Nsv,2*Nsv), State0(Nsv*2)
    DOUBLE PRECISION :: lhs(T), rhs(T,Nsv-1)

    ! NOTE: V0slopes (and thus also iV0slopes) is assumed to be block-diagonal, with separate blocks for each regression

    ! PART 1: Shockslopes
    offset = 0
    DO i=2,Nsv
       these = i - 1

       FORALL (k=1:T) lhs(k)       = y(i,k) / SVol(i,k)
       FORALL (k=1:T) rhs(k,1:i-1) = y(1:i-1,k) / SVol(i,k)

       ! slope indices are offset+1:offset+these
       call bayesRegressionSlope(shockslopes(offset+1:offset+these), lhs, rhs(:,1:i-1), these, T, 1.0d0, E0slopes(offset+1:offset+these), iV0slopes(offset+1:offset+these,offset+1:offset+these), VSLstream)

       FORALL (k=1:T) y(i,k) = lhs(k) * SVol(i,k)

       offset = offset + these

    END DO


    ! PART 2: SV

    ! log-linear observer
    logy2 = log(y ** 2 + 0.001d0)

    ! PART 2, STEP 1: DRAW KAI2STATES
    ! a) construct PDF for draws (stored in kai2CDF)
    FORALL (s=1:KSCmix,k=1:Nsv,j=1:T)
       kai2CDF(k,j,s) = exp(-0.5d0 * ((logy2(k,j) - h(k,j) - KSCmean(s)) / KSCvol(s))** 2) / KSCvol(s) * KSCpdf(s)
    END FORALL

    ! b) convert PDF into CDF for draws
    DO s=2,KSCmix
       kai2CDF(:,:,s) = kai2CDF(:,:,s-1) + kai2CDF(:,:,s)
    END DO
    DO s=1,KSCmix-1
       kai2CDF(:,:,s) = kai2CDF(:,:,s) / kai2CDF(:,:,KSCmix)
    END DO
    kai2CDF(:,:,KSCmix) = 1.0d0


    ! c) draw kai2states
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, Nsv * T, u, 0.0d0, 1.0d0 )
    FORALL (k=1:Nsv,j=1:T)
       kai2states(k,j) = COUNT(u(k,j) > kai2CDF(k,j,:)) + 1
    END FORALL


    ! PART 2, STEP 2: KALMAN FILTER FOR h

    ! prepare initial trend variance and noise variance
    FORALL (k=1:Nsv,j=1:T)
       volnoisemix(k,j) = KSCvol(kai2states(k,j))
    END FORALL

    ! demeaned observables 
    FORALL(k=1:Nsv,j=1:T)
       logy2star(k,j) = logy2(k,j) - KSCmean(kai2states(k,j))
    END FORALL

    ! state space matrices
    A = 0.0d0
    FORALL(k=1:Nsv,j=1:T) A(k,k,j) = rho
    FORALL(k=Nsv+1:Nsv*2,j=1:T) A(k,k,j) = 1.0d0

    B = 0.0d0
    FORALL(j=1:T) B(1:Nsv,1:Nsv,j) = sqrtVhshock

    C = 0.0d0
    FORALL(k=1:Nsv,j=1:T) C(k,k,j) = 1.0d0
    FORALL(k=1:Nsv,j=1:T) C(k,Nsv+k,j) = 1.0d0

    State0 = 0.0d0
    State0(Nsv+1:2*Nsv) = Eh0

    sqrtState0V = 0.0d0
    ! call eye(sqrtState0V, 10.0d0)
    sqrtState0V(Nsv+1:2*Nsv,Nsv+1:2*Nsv) = sqrtVh0

    CALL DLYAP(sqrtState0V(1:Nsv,1:Nsv), A(1:Nsv,1:Nsv,1), sqrtVhshock, Nsv, Nsv, errcode) 
    if (errcode /= 0) then
       write (*,*) 'DLYAP error (sqrtState0V)', errcode
       stop 1
    end if
    ! Factorize 
    CALL DPOTRF('L', Nsv, sqrtState0V(1:Nsv,1:Nsv), Nsv, errcode) ! recall: DLYAP returns fully symmetric matrix
    if (errcode /= 0) then
       write (*,*) 'DPOTRF error (sqrtState0V)', errcode
       stop 1
    end if
    ! zero out the upper triangular
    FORALL (i=1:Nsv-1) sqrtState0V(i,i+1:Nsv) = 0.0d0

    CALL samplerA3B3C3noise(State,StateShock,y2noise,logy2star,T,Nsv,Nsv*2,Nsv,A,B,C,volnoisemix,State0,sqrtState0V,VSLstream,errcode)
    if (errcode /= 0) then
       print *, 'something off with KSCvec sampler', errcode
       stop 1
    end if

    ! ! debug
    ! call savemat(A(:,:,1), 'A.debug')
    ! call savemat(B(:,:,1), 'B.debug')
    ! call savemat(C(:,:,1), 'C.debug')
    ! call savemat(State, 'State.debug')
    ! call savemat(StateShock, 'StateShock.debug')
    ! call savevec(State0, 'State0.debug')
    ! call savemat(sqrtState0V, 'sqrtState0V.debug')
    ! call savemat(logy2star, 'logy2star.debug')
    ! call savemat(volnoisemix, 'volnoisemix.debug')
    ! print *, 'rho', rho
    ! stop 33

    h      = State(1:Nsv,:) + State(Nsv+1:Nsv*2,:)
    hshock = StateShock(1:Nsv,:) 
    hbar   = State(Nsv+1:Nsv*2,0)

    SVol   = exp(h * 0.5d0)
  


  END SUBROUTINE SVHcholeskiKSCAR1corplus

  ! @\newpage\subsection{SVHdiffusecholeskiKSCAR1}@
  SUBROUTINE SVHdiffusecholeskiKSCAR1(T, Nsv, SVol, h, hshock, hbar, rho, shockslopes, Nslopes, y, sqrtVhshock, Eh0, sqrtVh0, VSLstream)

    ! h = log(SVol ** 2)

    INTENT(INOUT) :: SVol, VSLstream, y, h
    INTENT(IN) :: T,Nsv, sqrtVhshock, Eh0, sqrtVh0, rho
    INTENT(OUT) :: shockslopes, hshock, hbar

    DOUBLE PRECISION :: rho 

    INTEGER :: i, s, j, k, T, Nsv, Nslopes, offset, these, errcode
    DOUBLE PRECISION :: sqrtVhshock(Nsv,Nsv), Eh0(Nsv), sqrtVh0(Nsv,Nsv), shockslopes(Nslopes)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream

    ! KSC mixture
    INTEGER, PARAMETER :: KSCmix = 7
    DOUBLE PRECISION, DIMENSION(KSCmix), PARAMETER :: KSCmean = - 1.2704d0 + (/ -10.12999d0, -3.97281d0, -8.56686d0, 2.77786d0, .61942d0, 1.79518d0, -1.08819d0 /), KSCvar     = (/ 5.79596d0, 2.61369d0, 5.1795d0, .16735d0, .64009d0, .34023d0, 1.26261d0 /), KSCpdf= (/ .0073d0, .10556d0, .00002d0, .04395d0, .34001d0, .24566d0, .25750d0 /), KSCvol = sqrt(KSCvar)
    DOUBLE PRECISION, DIMENSION(Nsv,T,KSCmix) :: kai2CDF
    INTEGER, DIMENSION(Nsv,T) :: kai2states


    DOUBLE PRECISION, DIMENSION(Nsv,T) :: y, logy2, logy2star, volnoisemix, u
    DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: h, SVol
    DOUBLE PRECISION, DIMENSION(Nsv,1:T) :: hshock
    DOUBLE PRECISION, DIMENSION(Nsv) :: hbar

    ! state space matrices
    DOUBLE PRECISION :: A(Nsv*2,Nsv*2,T), B(Nsv*2,Nsv,T), C(Nsv,Nsv*2,T), y2noise(Nsv,T), State(Nsv*2,0:T), StateShock(Nsv*2,1:T), sqrtState0V(2*Nsv,2*Nsv), State0(Nsv*2)
    DOUBLE PRECISION :: lhs(T), rhs(T,Nsv-1)

    ! NOTE: V0slopes (and thus also iV0slopes) is assumed to be block-diagonal, with separate blocks for each regression

    ! PART 1: Shockslopes
    offset = 0
    DO i=2,Nsv
       these = i - 1

       FORALL (k=1:T) lhs(k)       = y(i,k) / SVol(i,k)
       FORALL (k=1:T) rhs(k,1:i-1) = y(1:i-1,k) / SVol(i,k)

       ! slope indices are offset+1:offset+these
       call bayesdiffuseRegressionSlope(shockslopes(offset+1:offset+these), lhs, rhs(:,1:i-1), these, T, 1.0d0, VSLstream)

       FORALL (k=1:T) y(i,k) = lhs(k) * SVol(i,k)

       offset = offset + these

    END DO


    ! PART 2: SV

    ! log-linear observer
    logy2 = log(y ** 2 + 0.001d0)

    ! PART 2, STEP 1: DRAW KAI2STATES
    ! a) construct PDF for draws (stored in kai2CDF)
    FORALL (s=1:KSCmix,k=1:Nsv,j=1:T)
       kai2CDF(k,j,s) = exp(-0.5d0 * ((logy2(k,j) - h(k,j) - KSCmean(s)) / KSCvol(s))** 2) / KSCvol(s) * KSCpdf(s)
    END FORALL

    ! b) convert PDF into CDF for draws
    DO s=2,KSCmix
       kai2CDF(:,:,s) = kai2CDF(:,:,s-1) + kai2CDF(:,:,s)
    END DO
    DO s=1,KSCmix-1
       kai2CDF(:,:,s) = kai2CDF(:,:,s) / kai2CDF(:,:,KSCmix)
    END DO
    kai2CDF(:,:,KSCmix) = 1.0d0


    ! c) draw kai2states
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, Nsv * T, u, 0.0d0, 1.0d0 )
    FORALL (k=1:Nsv,j=1:T)
       kai2states(k,j) = COUNT(u(k,j) > kai2CDF(k,j,:)) + 1
    END FORALL


    ! PART 2, STEP 2: KALMAN FILTER FOR h

    ! prepare initial trend variance and noise variance
    FORALL (k=1:Nsv,j=1:T)
       volnoisemix(k,j) = KSCvol(kai2states(k,j))
    END FORALL

    ! demeaned observables 
    FORALL(k=1:Nsv,j=1:T)
       logy2star(k,j) = logy2(k,j) - KSCmean(kai2states(k,j))
    END FORALL

    ! state space matrices
    A = 0.0d0
    FORALL(k=1:Nsv,j=1:T) A(k,k,j) = rho
    FORALL(k=Nsv+1:Nsv*2,j=1:T) A(k,k,j) = 1.0d0

    B = 0.0d0
    FORALL(j=1:T) B(1:Nsv,1:Nsv,j) = sqrtVhshock

    C = 0.0d0
    FORALL(k=1:Nsv,j=1:T) C(k,k,j) = 1.0d0
    FORALL(k=1:Nsv,j=1:T) C(k,Nsv+k,j) = 1.0d0

    State0 = 0.0d0
    State0(Nsv+1:2*Nsv) = Eh0

    sqrtState0V = 0.0d0
    ! call eye(sqrtState0V, 10.0d0)
    sqrtState0V(Nsv+1:2*Nsv,Nsv+1:2*Nsv) = sqrtVh0

    CALL DLYAP(sqrtState0V(1:Nsv,1:Nsv), A(1:Nsv,1:Nsv,1), sqrtVhshock, Nsv, Nsv, errcode) 
    if (errcode /= 0) then
       write (*,*) 'DLYAP error (sqrtState0V)', errcode
       stop 1
    end if
    ! Factorize 
    CALL DPOTRF('L', Nsv, sqrtState0V(1:Nsv,1:Nsv), Nsv, errcode) ! recall: DLYAP returns fully symmetric matrix
    if (errcode /= 0) then
       write (*,*) 'DPOTRF error (sqrtState0V)', errcode
       stop 1
    end if
    ! zero out the upper triangular
    FORALL (i=1:Nsv-1) sqrtState0V(i,i+1:Nsv) = 0.0d0

    CALL samplerA3B3C3noise(State,StateShock,y2noise,logy2star,T,Nsv,Nsv*2,Nsv,A,B,C,volnoisemix,State0,sqrtState0V,VSLstream,errcode)
    if (errcode /= 0) then
       print *, 'something off with KSCvec sampler', errcode
       stop 1
    end if

    ! ! debug
    ! call savemat(A(:,:,1), 'A.debug')
    ! call savemat(B(:,:,1), 'B.debug')
    ! call savemat(C(:,:,1), 'C.debug')
    ! call savemat(State, 'State.debug')
    ! call savemat(StateShock, 'StateShock.debug')
    ! call savevec(State0, 'State0.debug')
    ! call savemat(sqrtState0V, 'sqrtState0V.debug')
    ! call savemat(logy2star, 'logy2star.debug')
    ! call savemat(volnoisemix, 'volnoisemix.debug')
    ! print *, 'rho', rho
    ! stop 33

    h      = State(1:Nsv,:) + State(Nsv+1:Nsv*2,:)
    hshock = StateShock(1:Nsv,:) 
    hbar   = State(Nsv+1:Nsv*2,0)

    SVol   = exp(h * 0.5d0)
  

  END SUBROUTINE SVHdiffusecholeskiKSCAR1

! @\newpage\subsection{bayesRegressionSlope}@
  SUBROUTINE bayesRegressionSlope(bdraw, y, X, Nx, T, h, b0, V0i, VSLstream)
    ! draws vestor of regression slopes from Bayesian Regression 
    ! of scalar y on Vector X
    ! on exit, y return residuals

    INTENT(INOUT) :: y
    INTENT(IN) :: X, h, b0, V0i, Nx, T
    INTENT(OUT) :: bdraw
    INTENT(INOUT) :: VSLstream
    INTEGER :: T, Nx, status
    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: V0i, Vi
    DOUBLE PRECISION, DIMENSION(Nx) ::    b, b0, bdraw
    DOUBLE PRECISION :: h, y(T), X(T,Nx)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream


    ! prior/posterior

    ! b = inv(V0) * b0 + X'y * h
    call DSYMV('U', Nx, 1.0d0, V0i, Nx, b0, 1, 0.0d0, b, 1)
    call DGEMV('T', T, Nx, h, X, T, y, 1, 1.0d0, b, 1)

    ! Vi = inv(V0) + X'X * h
    Vi = V0i
    call DSYRK('U', 'T', Nx, T, h, X, T, 1.0d0, Vi, Nx)

    ! solve: Vi * b = inv(V0) * b0 + X'y * h
    call choleski(Vi) ! needed for inverting Vi as well as for draws, see below
    call DPOTRS('U', Nx, 1, Vi, Nx, b, Nx, status)
    if (status /= 0) then
       write(*,*) 'DPOTRS error: ', status, ' [BAYESREGRESSIONSLOPE]'
       stop 1
    end if
    
    ! draw from posterior with mean b and inverse variance Vi 
    ! notice: I am not scaling the draws by the choleski; using another factorization instead
    ! specifically: chol(Vi) * draws = z 
    ! (where V'= chol(Vl)' * chol(Vi), i.e. chol is upper triangular)
    status  = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx, bdraw, 0.0d0, 1.0d0)
    call DTRSV('U', 'N', 'N', Nx, Vi, Nx, bdraw, 1)

    ! add posterior mean
    bdraw = bdraw + b 
    
    ! resid = y - X * bdraw
    call dgemv('N', T, Nx, -1.0d0, X, T, bdraw, 1, 1.0d0, y, 1)

  END SUBROUTINE bayesRegressionSlope

! @\newpage\subsection{bayesDiffuseRegressionSlope}@
  SUBROUTINE bayesDiffuseRegressionSlope(bdraw, y, X, Nx, T, h, VSLstream)
    ! draws vestor of regression slopes from Bayesian Regression 
    ! of scalar y on Vector X
    ! on exit, y return residuals

    INTENT(INOUT) :: y
    INTENT(IN) :: X, h, Nx, T
    INTENT(OUT) :: bdraw
    INTENT(INOUT) :: VSLstream
    INTEGER :: T, Nx, status
    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: Vi
    DOUBLE PRECISION, DIMENSION(Nx) ::    b, bdraw
    DOUBLE PRECISION :: h, y(T), X(T,Nx)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream


    ! prior/posterior

    ! b = X'y * h
    ! call DSYMV('U', Nx, 1.0d0, V0i, Nx, b0, 1, 0.0d0, b, 1)
    b = 0.0d0
    call DGEMV('T', T, Nx, h, X, T, y, 1, 0.0d0, b, 1)

    ! Vi = X'X * h
    Vi = 0.0d0
    call DSYRK('U', 'T', Nx, T, h, X, T, 0.0d0, Vi, Nx)

    ! solve: Vi * b = inv(V0) * b0 + X'y * h
    call choleski(Vi) ! needed for inverting Vi as well as for draws, see below
    call DPOTRS('U', Nx, 1, Vi, Nx, b, Nx, status)
    if (status /= 0) then
       write(*,*) 'DPOTRS error: ', status, ' [BAYESREGRESSIONSLOPE]'
       stop 1
    end if
    
    ! draw from posterior with mean b and inverse variance Vi 
    ! notice: I am not scaling the draws by the choleski; using another factorization instead
    ! specifically: chol(Vi) * draws = z 
    ! (where V'= chol(Vl)' * chol(Vi), i.e. chol is upper triangular)
    status  = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx, bdraw, 0.0d0, 1.0d0)
    call DTRSV('U', 'N', 'N', Nx, Vi, Nx, bdraw, 1)

    ! add posterior mean
    bdraw = bdraw + b 
    
    ! resid = y - X * bdraw
    call dgemv('N', T, Nx, -1.0d0, X, T, bdraw, 1, 1.0d0, y, 1)

  END SUBROUTINE bayesDiffuseRegressionSlope

! @\newpage\subsection{bayesVARSV}@
  SUBROUTINE bayesVARSV(bdraw, Y, p, Ydata, Ny, T, iSigmaResid, b0, V0i, VSLstream)
    ! Bayesian VAR with (known) Time-varying Volatility
    ! note: iSigmaResid is (Ny,Ny,T) and assumed upper triangular
    ! (assuming no constant)
    ! on exit, Y returns residuals Y - X * reshape(bdraw, Nx, Ny)
    ! notice: top rows of companion have transpose(reshape(bdraw,Nx,Ny)) 

    INTENT(IN) :: Ydata, Ny, p, T, iSigmaResid, b0, V0i
    INTENT(OUT) :: Y, bdraw
    INTENT(INOUT) :: VSLstream
    INTEGER :: T, Nx, Ny, status, Nb, j, p
    DOUBLE PRECISION, DIMENSION(Ny * Ny * p,Ny * Ny * p) :: V0i, Vi
    DOUBLE PRECISION, DIMENSION(Ny * Ny * p) ::    b, b0, bdraw
    DOUBLE PRECISION :: iSigmaResid(Ny,Ny,T), Y(T,Ny), Ytilde(T,Ny), Ydata(-(p-1):T,Ny), X(T,Ny * p), XX(Ny * p, Ny * p), XY(Ny * p, Ny)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream
    
    ! STEP 1: construct regressors


    Y  = Ydata(1:T,:)
    Nx = Ny * p
    Nb = Ny * Nx
    FORALL (j = 1:p) X(:, ((j-1) * Ny + 1) : ((j-1) * Ny + Ny) ) = Ydata(-(j-1):T-j,:)

    ! STEP 2: estimate VAR coefficients: beta = inv(XX) * Xy
    ! construct Ytilde(t) = iSigmaResid(t) Y(t)
    Ytilde = 0.0d0
    DO j = 1, T
       call dsymv('U', Ny, 1.0d0, iSigmaResid(:,:,j), Ny, Y(j,:), 1, 0.0d0, Ytilde(j,:), 1)
    END DO
    ! call savemat(Y, 'y.debug')
    ! call savemat(Ytilde, 'ytilde.debug')
    ! call savemat(X, 'x.debug')
    
    ! XY = sum_t {X(t) Ytilde(t)'}
    call DGEMM('T', 'N', Nx, Ny, T, 1.0d0, X, T, Ytilde, T, 0.0d0, XY, Nx)
    ! store vec(XY) in b
    call vec(b,XY)

    ! b = V0i * b0 + b 
    ! (b is not yet complete, need to multiply by posterior Variance, see below)
    call DSYMV('U', Nb, 1.0d0, V0i, Nb, b0, 1, 1.0d0, b, 1)

    ! POSTERIOR VARIANCE
    ! Vi = V0i + sum_t kron(iSigmaResid(t), X(t) X(t)')

    Vi = V0i
    DO j = 1, T

       XX = 0.0d0
       call DSYR('U',Nx,1.0d0,X(j,:),1,XX,Nx)
       call symmetric(XX) ! important for symkronecker
       call symkronecker(1.0d0,iSigmaResid(:,:,j),Ny,XX,Nx,1.0d0,Vi)
    END DO
    


    
    ! Solve for posterior mean
    ! solve: Vi * b = ...
    ! 1) Choleski factorization
    call dpotrf('u', Nb, Vi, Nb, status)
    if (status /= 0) then
       write(*,*) 'CHOLESKI ERROR:', status, ' [BAYESVAR]'

       ! do j=1,T

       ! end do
       

       call savemat(V0i, 'V0i.dat.debug') ! debug
       call savemat(Vi, 'Vi.dat.debug') ! debug
       call savevec(b0, 'b0.dat.debug') ! debug
       call savevec(b, 'b.dat.debug') ! debug
       call savemat(X, 'X.dat.debug') ! debug
       call savemat(Y, 'Y.dat.debug') ! debug


       stop 1
    end if

    ! zero out lower triangular
    forall (j = 1 : Nb-1) Vi(j+1:Nb,j) = 0.0d0
    
    ! 2) solve Vi * b = z
    call DPOTRS('U', Nb, 1, Vi, Nb, b, Nb, status)
    if (status /= 0) then
       write(*,*) 'DPOTRS error: ', status, ' [BAYESVAR]'
       stop 1
    end if

    ! DRAW FROM POSTERIOR with mean b and inverse variance Vi 
    ! notice: I am not scaling the draw)' * chol(Vi), i.e. chol is upper triangular)
    status  = vdrnggaussian(VSLmethodGaussian, VSLstream, Nb, bdraw, 0.0d0, 1.0d0)
    call DTRSV('U', 'N', 'N', Nb, Vi, Nb, bdraw, 1)

    ! add posterior mean
    bdraw = bdraw + b
    
    ! resid = Y - X * reshape(beta, Nx, Ny)
    ! note: DGEMM does not care about explicitly reshaping beta
    call DGEMM('N','N',T,Ny,Nx,-1.0d0,X,T,bdraw,Nx,1.0d0,Y,T)

    
  END SUBROUTINE bayesVARSV

! @\newpage\subsection{bayesdiffuseVARSV}@
  SUBROUTINE bayesdiffuseVARSV(bdraw, Y, p, Ydata, Ny, T, iSigmaResid, VSLstream)
    ! Bayesian VAR with (known) Time-varying Volatility
    ! note: iSigmaResid is (Ny,Ny,T) and assumed upper triangular
    ! (assuming no constant)
    ! on exit, Y returns residuals Y - X * reshape(bdraw, Nx, Ny)
    ! notice: top rows of companion have transpose(reshape(bdraw,Nx,Ny)) 

    INTENT(IN) :: Ydata, Ny, p, T, iSigmaResid
    INTENT(OUT) :: Y, bdraw
    INTENT(INOUT) :: VSLstream
    INTEGER :: T, Nx, Ny, status, Nb, j, p
    DOUBLE PRECISION, DIMENSION(Ny * Ny * p,Ny * Ny * p) :: Vi
    DOUBLE PRECISION, DIMENSION(Ny * Ny * p) ::    b, bdraw
    DOUBLE PRECISION :: iSigmaResid(Ny,Ny,T), Y(T,Ny), Ytilde(T,Ny), Ydata(-(p-1):T,Ny), X(T,Ny * p), XX(Ny * p, Ny * p), XY(Ny * p, Ny)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    TYPE (vsl_stream_state) :: VSLstream
    
    ! STEP 1: construct regressors


    Y  = Ydata(1:T,:)
    Nx = Ny * p
    Nb = Ny * Nx
    FORALL (j = 1:p) X(:, ((j-1) * Ny + 1) : ((j-1) * Ny + Ny) ) = Ydata(-(j-1):T-j,:)

    ! STEP 2: estimate VAR coefficients: beta = inv(XX) * Xy
    ! construct Ytilde(t) = iSigmaResid(t) Y(t)
    Ytilde = 0.0d0
    DO j = 1, T
       call dsymv('U', Ny, 1.0d0, iSigmaResid(:,:,j), Ny, Y(j,:), 1, 0.0d0, Ytilde(j,:), 1)
    END DO
    ! call savemat(Y, 'y.debug')
    ! call savemat(Ytilde, 'ytilde.debug')
    ! call savemat(X, 'x.debug')
    
    ! XY = sum_t {X(t) Ytilde(t)'}
    call DGEMM('T', 'N', Nx, Ny, T, 1.0d0, X, T, Ytilde, T, 0.0d0, XY, Nx)
    ! store vec(XY) in b
    call vec(b,XY)

    ! b = V0i * b0 + b 
    ! (b is not yet complete, need to multiply by posterior Variance, see below)
    ! call DSYMV('U', Nb, 1.0d0, V0i, Nb, b0, 1, 1.0d0, b, 1)
   
    ! POSTERIOR VARIANCE
    ! Vi = V0i + sum_t kron(iSigmaResid(t), X(t) X(t)')

    Vi = 0.0d0
    DO j = 1, T

       XX = 0.0d0
       call DSYR('U',Nx,1.0d0,X(j,:),1,XX,Nx)
       call symmetric(XX) ! important for symkronecker
       call symkronecker(1.0d0,iSigmaResid(:,:,j),Ny,XX,Nx,1.0d0,Vi)
    END DO
    


    
    ! Solve for posterior mean
    ! solve: Vi * b = ...
    ! 1) Choleski factorization
    call dpotrf('u', Nb, Vi, Nb, status)
    if (status /= 0) then
       write(*,*) 'CHOLESKI ERROR:', status, ' [BAYESVAR]'

       call savemat(Vi, 'Vi.dat.debug') ! debug
       call savevec(b, 'b.dat.debug') ! debug
       call savemat(X, 'X.dat.debug') ! debug
       call savemat(Y, 'Y.dat.debug') ! debug
       call savearray3(iSigmaResid, 'iSigmaResid', 'debug') 

       stop 1
    end if

    ! zero out lower triangular
    forall (j = 1 : Nb-1) Vi(j+1:Nb,j) = 0.0d0
    
    ! 2) solve Vi * b = z
    call DPOTRS('U', Nb, 1, Vi, Nb, b, Nb, status)
    if (status /= 0) then
       write(*,*) 'DPOTRS error: ', status, ' [BAYESVAR]'
       stop 1
    end if

    ! DRAW FROM POSTERIOR with mean b and inverse variance Vi 
    ! notice: I am not scaling the draw)' * chol(Vi), i.e. chol is upper triangular)
    status  = vdrnggaussian(VSLmethodGaussian, VSLstream, Nb, bdraw, 0.0d0, 1.0d0)
    call DTRSV('U', 'N', 'N', Nb, Vi, Nb, bdraw, 1)

    ! add posterior mean
    bdraw = bdraw + b
    
    ! resid = Y - X * reshape(beta, Nx, Ny)
    ! note: DGEMM does not care about explicitly reshaping beta
    call DGEMM('N','N',T,Ny,Nx,-1.0d0,X,T,bdraw,Nx,1.0d0,Y,T)

    
  END SUBROUTINE bayesdiffuseVARSV
! @\newpage\subsection{VARmaxroot}@
  SUBROUTINE VARmaxroot(maxlambda, beta, Ny, p)
    
    ! assumes beta contains no slope for constant

    INTENT(IN) :: beta, Ny, p
    INTENT(OUT) :: maxlambda

    INTEGER :: Ny, p, Nx, j
    DOUBLE PRECISION  :: kompanion(Ny * p, Ny * p), beta(Ny * p, Ny), maxlambda


    ! construct companion form
    Nx = Ny * p
    kompanion = 0.0d0
    FORALL (j = 1 : Ny * (p - 1))
       kompanion(j+Ny,j) = 1.0d0
    END FORALL

    kompanion(1:Ny,:) = transpose(beta(1:Nx,:))

    ! compute largest root


    maxlambda = maxroot(kompanion, Nx)

  END SUBROUTINE VARmaxroot

! @\newpage\subsection{ARmaxroot}@
  SUBROUTINE ARmaxroot(maxlambda, beta, p)
    
    ! assumes beta contains no slope for constant

    INTENT(IN) :: beta, p
    INTENT(OUT) :: maxlambda

    INTEGER :: p, j
    DOUBLE PRECISION  :: kompanion(p, p), beta(p), maxlambda


    ! construct companion form
    kompanion = 0.0d0
    FORALL (j = 1 : (p - 1))
       kompanion(1+j,j) = 1.0d0
    END FORALL

    kompanion(1,:) = beta

    ! compute largest root


    maxlambda = maxroot(kompanion, p)

  END SUBROUTINE ARmaxroot
 
! @\newpage\subsection{simPriorMaxroot}@
  SUBROUTINE simPriorMaxroot(maxlambdas, Ndraws, f0, sqrtVf0, Ny, p, VSLstream)
    
    ! assumes beta contains slope on constant

    INTENT(IN) :: f0, sqrtVf0, Ny, p, Ndraws
    INTENT(OUT) :: maxlambdas
    INTENT(INOUT) :: VSLstream

    INTEGER :: Ny, p, Ndraws, Nx, errcode, j, Nf
    DOUBLE PRECISION  :: maxlambdas(Ndraws), kompanion(Ny * p, Ny * p), f(Ny * Ny * p,Ndraws), f0(Ny * Ny * p), sqrtVf0(Ny * Ny * p, Ny * Ny * p)
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    type (vsl_stream_state) :: VSLstream


    Nf = Ny * Ny * p

    ! construct companion form
    Nx = Ny * p
    kompanion = 0.0d0
    FORALL (j = 1 : Ny * (p - 1))
       kompanion(j+Ny,j) = 1.0d0
    END FORALL

    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Ndraws * Nf, f, 0.0d0, 1.0d0)

    DO j=1,Ndraws

       CALL DTRMV('U', 'T', 'N', Nf, sqrtVf0, Nf, f(:,j), 1)
       f(:,j) = f(:,j) + f0

       kompanion(1:Ny,:) = transpose(reshape(f(:,j), (/ Ny * p, Ny /)))
       
       maxlambdas(j) = maxroot(kompanion, Nx)

    END DO


  END SUBROUTINE simPriorMaxroot
END MODULE gibbsbox


