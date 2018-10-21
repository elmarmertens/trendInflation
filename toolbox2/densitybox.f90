MODULE densitybox

  USE embox, only : sqrttwopi, logtwopi, mean, variance, normpdf, savemat, savevec, savearray3
  USE blaspack, only: sandwich !, sandwich2
  USE statespacebox, only : simyABCsvol0

  USE vslbox


  IMPLICIT NONE

CONTAINS

! @\newpage\subsection{predictiveDensitySV}@
SUBROUTINE predictiveDensitySV(ypdf, yforecast, ycondvar, Nhorizons, horizons, maxhorizons, Ny, ypred, ynanpred, Nx, Nw, x, A, B, C, Nsv, h0, hbar, hrho, hinno, VSLstream)


  INTENT(INOUT) :: VSLstream
  INTENT(OUT)   :: ypdf, yforecast, ycondvar
  INTENT(IN)    :: Nhorizons, horizons, maxhorizons, Ny, ypred, ynanpred, Nx, Nw, x, A, B, C, Nsv, h0, hbar, hrho, hinno

  INTEGER, PARAMETER :: NsimSV = 100
  INTEGER :: Nhorizons, horizons(Nhorizons), maxhorizons
  INTEGER :: Ny, Nx, Nw, Nsv

  double precision :: x(Nx), A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)
  double precision, dimension(Nsv) :: h0, hbar, hrho, hinno

  double precision, dimension(Ny, Nhorizons) :: ypdf, yforecast, ycondvar
  double precision, dimension(Ny, maxhorizons) :: ypred
  logical, dimension(Ny, maxhorizons) :: ynanpred

  DOUBLE PRECISION, dimension(Ny,maxhorizons) :: forecastMean  
  DOUBLE PRECISION, dimension(NsimSV,Ny,maxhorizons) :: forecastCondVar

  INTEGER :: hh, ii, jj, kk, nn

  ! variables for generating mean forecasts
  double precision, dimension(Nx) :: xlag, xhat


  ! helper
  DOUBLE PRECISION, DIMENSION(Nsv,maxhorizons,NsimSV) :: h, hshock
  DOUBLE PRECISION :: SigmaY(Ny,Ny,maxhorizons), VMA(Nx,Nw,0:maxhorizons-1)
  
  type (vsl_stream_state) :: VSLstream
  integer :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


  ! init
  ypdf             = 0.0d0
  yforecast        = 0.0d0
  ycondvar         = 0.0d0

  forecastMean     = 0.0d0
  forecastCondVar  = 0.0d0

  ! FORECAST MEAN
  xhat=x ! necessary copy since x is intent(in)

  DO hh = 1, maxhorizons

     xlag = xhat
     CALL DGEMV('N', Nx, Nx, 1.0d0, A,  Nx, xlag, 1, 0.0d0, xhat, 1)
     CALL DGEMV('N', Ny, Nx, 1.0d0, C,  Ny, xhat, 1, 0.0d0, forecastMean(:,hh), 1)

  END DO ! hh 


  ! TODO: rewrite a separate RW-SV version (code below works also for rho=1)

  ! SIMULATE SVs (needed for variance computations) 
  ! draw random numbers
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * maxhorizons * NsimSV, hshock, 0.0d0, 1.0d0)
  ! compute transition from jump off (demeaned)
  hh = 1
  FORALL (jj=1:Nsv,ii=1:NsimSV) h(jj,hh,ii)    = hrho(jj) * h0(jj)        + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii)
  DO hh = 2,maxhorizons
     FORALL (jj=1:Nsv,ii=1:NsimSV) h(jj,hh,ii) = hrho(jj) * h(jj,hh-1,ii) + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii) ! for these loops, hh should ideally be in third dimension, but next set of loops works better with ii in last dimension
  END DO

  VMA        = 0.0d0
  VMA(:,:,0) = B
  DO hh = 1,maxhorizons-1
      call dgemm('N', 'N', Nx, Nw, Nx, 1.0d0, A, Nx, VMA(:,:,hh-1), Nx, 0.0d0, VMA(:,:,hh), Nx)
  END DO


  ! FORECAST VARIANCES (unconstrained) FOR EACH SIMULATED SV
  DO ii =1,NsimSV
     ! SigmaY = condvarABCSV(A,B,C,h(:,:,ii),Ny,Nx,Nw,Nsv,maxhorizons) 
     SigmaY = condvarVMACSV(VMA,C,h(:,:,ii),Ny,Nx,Nw,Nsv,maxhorizons) 
     forall (nn=1:Ny,hh=1:maxhorizons) forecastCondVar(ii,nn,hh) = SigmaY(nn,nn,hh) ! not quite efficient to have ii as the first-column argument in *this* loop, but helps with mean() further below
  END DO


  ! call savemat(Sigma, 'Sigma.debug')
  ! call savemat(SigmaY, 'SigmaY.debug')
  ! call savemat(B, 'B.debug')
  ! call savemat(A, 'A.debug')
  ! call savemat(C, 'C.debug')
  ! call savemat(forecastcondvar(1,:,:), 'condvar.debug')
  ! call savearray3(VMA, 'VMA', 'debug')
  ! call savemat(h(:,:,1), 'h.debug')
  ! stop 1

  ! yforecast
  forall (nn=1:Ny,kk=1:Nhorizons) yforecast(nn,kk) = forecastMean(nn,horizons(kk))
  ! ycondvar
  forall (nn=1:Ny,kk=1:Nhorizons) ycondvar(nn,kk)  = mean(forecastCondVar(:,nn,horizons(kk))) ! Note: there is no variance across conditional expectations here

  ! compute llf 
  do kk=1,Nhorizons
     hh = horizons(kk)
     do nn = 1, Ny ! could replace this with a where constuct, but might make code less obvious
        ! pdf (if predicted observations available)
        if (.not. ynanpred(nn,hh)) then
           ! normal pdf
           ypdf(nn,kk) = mean(exp(-0.5d0 * (((ypred(nn,hh) - forecastMean(nn,hh)) ** 2) / forecastCondVar(:,nn,hh))) /  (sqrttwopi * sqrt(forecastCondVar(:,nn,hh))))
        end if

     end do ! nn 
  end do ! kk


END SUBROUTINE predictiveDensitySV

! @\newpage\subsection{condvarVMACSV}@
FUNCTION condvarVMACSV(VMA,C,hSV,Ny,Nx,Nw,Nsv,horizons) result(SigmaY)

  INTENT(IN) :: VMA,C,hSV,Ny,Nx,Nw,Nsv,horizons

  integer :: Ny,Nx,Nw,Nsv,horizons
  double precision :: VMA(Nx,Nw,0:horizons-1), C(Ny,Nx), hSV(Nsv,horizons)
  double precision :: SigmaY(Ny,Ny,horizons)
  double precision :: VMAsv(Nx,Nw), SigmaX(Nx,Nx)

  integer jj,hh,ll,pit

  SigmaY = 0.0d0
  DO hh=1,horizons
     SigmaX  = 0.0d0

     DO ll=0,hh-1
        ! rescale VMA to account for SV components
        VMAsv  = VMA(:,:,ll)
        pit    = hh-ll 
        ! scale SV columns
        FORALL (jj=1:Nsv) VMAsv(:,jj) = VMAsv(:,jj) * exp(0.5d0 * hSV(jj,pit)) 
        ! Sigma = VMAsv * VMAsv' + Sigma
        call dsyrk('U', 'N', Nx, Nw, 1.0d0, VMAsv, Nx, 1.0d0, SigmaX, Nx)
     END DO
     call sandwich(SigmaY(:,:,hh), C, Ny, SigmaX, Nx)

  END DO


END FUNCTION condvarVMACSV

END MODULE densitybox


