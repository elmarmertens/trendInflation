MODULE otherbox

  IMPLICIT NONE

CONTAINS

!// Kronecker.f90:
!//
!// Given matrices A and B, the FUNCTION 'Kronecker' in below calculates
!//				C = alpha*C + beta*(A kron B)
!// It is generic Fortran90 code that does not require any installation of other libraries.
!// One is able to use it on the fly.
!// For simple instructions please read comments.
!//
!// author: Feng Chen, chenfeng3372338@gmail.com.
!// date: 2012/01/06.
!// version: 2.8.
!//

subroutine Kronecker(ma, na, mb, nb, mc, nc, alpha, beta, A, B, C)
	integer, intent(in) :: ma, na, mb, nb
	real(kind=8), intent(in) :: alpha, beta
	real(kind=8), dimension(ma,na), intent(in) :: A
	real(kind=8), dimension(mb,nb), intent(in) :: B
	real(kind=8), dimension(mc,nc), intent(inout) :: C
	
	integer :: i, j, i1, j1  
	
	
	do j = 1, na
		do i = 1, ma
			i1 = (i-1)*mb + 1
			j1 = (j-1)*nb + 1 
			C(i1:i1+mb-1,j1:j1+nb-1) = alpha*C(i1:i1+mb-1,j1:j1+nb-1) + beta*A(i,j)*B 
		end do
	end do
			
end subroutine Kronecker 



END MODULE otherbox
