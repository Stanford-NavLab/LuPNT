c
c	rotation matrix t3=<-psi,X>
c
c	iverse = 1 straight transform
c	      =-1 inverse transform
c
	subroutine t3(itime,x_in,x_out,iverse)
	real x_in(3),x_out(3),q_c(3),psi
	integer*4 itime(2),iverse
c
	call get_q_c(itime,q_c)
c
	if(q_c(3).eq.0.0) then
	  psi=-sign(1.5707963,q_c(2))
	else
	  psi=-atan(q_c(2)/q_c(3))
	end if
	call rotate_x(iverse*psi,x_in,x_out)
c
	return
	end
