subroutine elemsubr ( npelm, x, y, nunk_pel, elem_mat,  &
                      elem_vec, elem_mass, prevsolution, itype )

!                  INPUT / OUTPUT PARAMETERS

implicit none
integer, intent(in) :: npelm, nunk_pel, itype
double precision, intent(in) :: x(1:npelm), y(1:npelm),  &
                                prevsolution(1:nunk_pel)
double precision, intent(out) :: elem_mat(1:nunk_pel,1:nunk_pel),  &
                                 elem_vec(1:nunk_pel),  &
                                 elem_mass(1:nunk_pel)

!******************************************************************
!
!                  LOCAL PARAMETERS
!
!==================================================================

double precision :: beta(1:3), gamma(1:3), delta, h 
double precision :: pi = 4d0*atan(1d0), rho=1d0, alpha = 100d0 
integer :: i, j

if (itype==1) then

  delta = (x(2)-x(1))*(y(3)-y(1)) - (y(2)-y(1))*(x(3)-x(1))
  beta(1) = (y(2)-y(3))/delta
  beta(2) = (y(3)-y(1))/delta
  beta(3) = (y(1)-y(2))/delta
  gamma(1) = (x(3)-x(2))/delta
  gamma(2) = (x(1)-x(3))/delta
  gamma(3) = (x(2)-x(1))/delta

  do j = 1,3
    do i = 1,3
      elem_mat(i,j)=(pi/3)*abs(delta)* &
        (beta(i)*beta(j) + gamma(i)*gamma(j))*sum(x)
    enddo
  enddo

  elem_vec(1:3) = 0d0

elseif(itype==2) then ! upper boundary

  h = sqrt((x(2)-x(1))**2 + (y(2)-y(1))**2)

  do j = 1,2
    do i = 1,2
      if(i==j) then
        elem_mat(i,j)=pi*alpha*h*x(i)
      else 
        elem_mat(i,j)=0d0
      endif
    enddo
  enddo
  
  elem_vec(1:2) = pi*h*x*y/rho

elseif(itype==3) then  ! lower boundary 

  h = sqrt((x(2)-x(1))**2 + (y(2)-y(1))**2)

  do j = 1,2
    do i = 1,2
      if(i==j) then
        elem_mat(i,j)=pi*alpha*h*x(i)
      else 
        elem_mat(i,j)=0d0
      endif
    enddo
  enddo

  elem_vec(1:2) = 0d0
  
endif
end subroutine elemsubr
