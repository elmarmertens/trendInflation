! include 'mkl_vsl.fi'

MODULE vslbox

! USE mkl_vsl
! USE mkl_vsl_type

IMPLICIT NONE

type vsl_stream_state
   integer :: state
end type vsl_stream_state

integer, parameter :: vsl_brng_mt19937=0, vsl_brng_mt2203=0
CONTAINS

INTEGER FUNCTION vslnewstream(VSLstream, brng, seed)     
INTENT (INOUT) :: VSLstream
INTENT (IN) :: brng, seed
INTEGER :: brng, seed, n
type (vsl_stream_state) :: VSLstream
INTEGER, DIMENSION(:), ALLOCATABLE :: gseed

CALL RANDOM_SEED(size = n)
ALLOCATE(gseed(n))
gseed = seed
CALL RANDOM_SEED(PUT = gseed)
DEALLOCATE(gseed)

VSLstream%state = 0 
vslnewstream = 0


END FUNCTION vslnewstream

INTEGER FUNCTION vsldeletestream(VSLstream)     
INTENT (INOUT) :: VSLstream
type (vsl_stream_state) :: VSLstream

vsldeletestream = 0
END FUNCTION vsldeletestream

INTEGER FUNCTION vdrnggaussian(VSLmethod, VSLstream, T, x, mu, sig)     
INTENT (INOUT) :: VSLstream, x
INTENT (IN) :: VSlmethod, T, mu, sig
INTEGER :: VSLmethod, T
DOUBLE PRECISION :: mu, sig
type (vsl_stream_state) :: VSLstream
DOUBLE PRECISION, DIMENSION(T) :: x
DOUBLE PRECISION, DIMENSION(T,12) :: u

vdrnggaussian = 0

CALL RANDOM_NUMBER(u)
x = sum(u,2) - 6

x = x * sig + mu

END FUNCTION vdrnggaussian


INTEGER FUNCTION vdrnguniform(VSLmethod, VSLstream, T, u, a, b)     
INTENT (INOUT) :: VSLstream, u
INTENT (IN) :: VSlmethod, T, a, b
INTEGER :: VSLmethod, T
DOUBLE PRECISION :: a, b
type (vsl_stream_state) :: VSLstream
DOUBLE PRECISION, DIMENSION(T) :: u

vdrnguniform = 0

CALL RANDOM_NUMBER(u)
u = a + (b - a) * u

END FUNCTION vdrnguniform




END MODULE vslbox
