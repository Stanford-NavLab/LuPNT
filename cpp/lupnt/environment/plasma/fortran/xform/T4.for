c
c	transformation matrix t4=<-mu,Y>
c
c	iverse = 1 forward transform
c	       =-1 inverse transform
c
	subroutine t4 (itime,x_in,x_out,iverse)
	real x_in(3),x_out(3),mu,q_c(3)
	integer*4 itime(2),iverse
c
	call get_q_c(itime,q_c)
c
	mu=-atan(q_c(1)/sqrt(q_c(2)*q_c(2)+q_c(3)*q_c(3)))
	call rotate_y(iverse*mu,x_in,x_out)
c
	return
	end
