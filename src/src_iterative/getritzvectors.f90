! Copyright (C) 2015-2023 exciting team (Berlin and Riga)
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getritzvectors(l,m,n,A,B,C,nstart,npw)
!>  Constructs the trial wavefunctions after subspace diagonalisation. 
!> 
      Use constants, Only : zzero, zone
      use mod_gkvector
      use modmpi
      implicit none
      integer, intent(in) :: l           ! leading dimension, size of the basis       
      integer, intent(in) :: m           ! size of the subspace 
      integer, intent(in) :: n           ! number of Ritz eigenpairs 
      integer, intent(in) :: nstart      ! offset in the subspace; typically used to skip the subspace vectors that correspond strictly to LOs
      integer, intent(in) :: npw         ! number of APWs (not used in the default version)
      complex(8), intent(in) :: A(l,m)   ! subspace vectors
      complex(8), intent(in) :: B(m,m)   ! matrix of eigenvectors
      complex(8), intent(out) :: C(l,n)  ! new trial wavefunctions

      real(8) :: ta,tb
      integer :: nlo,i,j
      integer, allocatable :: offset(:),ngklist(:),ibuf(:)

call timesec(ta)
if (.true.) then
! Default option:
! we keep it simple and calculate the trial wavefunctions in a single go
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  l, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  m, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  A, &           ! A
                  l,&           ! LDA ... leading dimension of A
                  B(1,nstart), &           ! B
                  m, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  C, &  ! C
                  l &      ! LDC ... leading dimension of C
                  )
else  
! Disabled option (not tested for a while):
! calculate the APW and LO parts separately 
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  npw, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  m, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  A(1,1), &           ! A
                  l,&           ! LDA ... leading dimension of A
                  B(1,nstart), &           ! B
                  m, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  C(1,1), &  ! C
                  l &      ! LDC ... leading dimension of C
                  )
      nlo=l-npw
      if (nlo.ne.0) then
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    nlo, &          ! M ... rows of op( A ) = rows of C
                    n, &           ! N ... cols of op( B ) = cols of C
                    m, &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    A(npw+1,1), &           ! A
                    l,&           ! LDA ... leading dimension of A
                    B(1,nstart), &           ! B
                    m, &          ! LDB ... leading dimension of B
                    zzero, &          ! beta
                    C(npw+1,1), &  ! C
                    l &      ! LDC ... leading dimension of C
                    )
     endif

endif


!endif

call timesec(tb)
#ifdef TIMINGS
write(*,*) 'getritzvectors', tb-ta
#endif
end subroutine getritzvectors
