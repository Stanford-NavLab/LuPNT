c
c	transformation matrix t1=<theta,z>
c
c	iverse = 1 straight transform
c	      =-1 inverse transform
c
	subroutine t1(itime,x_in,x_out,iverse)
c
	real theta,x_in(3),x_out(3),t0,ut,temp
	integer*4 itime(2),iday,iyr,iverse
	parameter (degrad=3.1415927/180.0)
c
	temp=t0(itime,iyr,iday,ut)
	theta=(100.461+36000.770*temp+15.04107*ut)*degrad
c
	call rotate_z (float(iverse)*theta,x_in,x_out)
c
	return
	end
