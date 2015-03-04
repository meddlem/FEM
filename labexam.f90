      program labexam

!     --- Sample program for the lab numerical methods in scientific computing

      call startsepcomp
      
      end program labexam
      
      subroutine elemsubr ( npelm, x, y, nunk_pel, elem_mat,  &
                            elem_vec, elem_mass, prevsolution, itype )
!
!                       INPUT / OUTPUT PARAMETERS
!
      implicit none
      integer, intent(in) :: npelm, nunk_pel, itype
      double precision, intent(in) :: x(1:npelm), y(1:npelm),  &
                                      prevsolution(1:nunk_pel)
      double precision, intent(out) :: elem_mat(1:nunk_pel,1:nunk_pel),  &
                                       elem_vec(1:nunk_pel),  &
                                       elem_mass(1:nunk_pel)
      
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision :: mu, beta(1:3), gamma(1:3), delta
      integer :: i, j

! ======================================================================
!
      if ( itype == 1 ) then

!     --- Type = 1: mu = 1

         mu = 1d0

      else

!     --- Type = 2: mu = 2

         mu = 2d0

      end if

!     --- Compute the factors e_ij and delta as defined in the manual

      delta = (x(2)-x(1))*(y(3)-y(1))-(y(2)-y(1))*(x(3)-x(1))
      beta(1) = (y(2)-y(3))/delta
      beta(2) = (y(3)-y(1))/delta
      beta(3) = (y(1)-y(2))/delta
      gamma(1) = (x(3)-x(2))/delta
      gamma(2) = (x(1)-x(3))/delta
      gamma(3) = (x(2)-x(1))/delta

!     --- Fill the element matrix as defined in the Lecture Notes

      do j = 1, 3
         do i = 1, 3
            elem_mat(i,j) = mu * 0.5d0 * abs(delta) *   &
               ( beta(i)*beta(j) + gamma(i)*gamma(j) )
         end do
      end do

!     --- The element vector is zero

      elem_vec(1:3) = 0d0

      end subroutine elemsubr
