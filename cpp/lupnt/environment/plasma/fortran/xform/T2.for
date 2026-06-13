c
c	transformation matrix t2=<lamdas,Z>*<epsilon,X>
c
c	iverse = 1 straight transform
c	      =-1 inverse transform
c
	subroutine t2(itime,x_in,x_out,iverse)
	real x_in(3),x_out(3),m,cgamma,lamdas,epsilon
	real temp(3),ut,tt0,t0
	integer*4 itime(2),iyr,iday,iverse
	parameter (degrad=3.1415927/180.0)
c
	tt0=t0(itime,iyr,iday,ut)
	epsilon=(23.439-0.013*tt0)*degrad
	m=(357.528+35999.05*tt0+0.04107*ut)*degrad
	cgamma=280.46+36000.772*tt0+0.04107*ut
	lamdas=(cgamma+(1.915-0.0048*tt0)*sin(m)+0.02*sin(2.0*m))*degrad
c
	if (iverse.eq.1) then
	  call rotate_x (epsilon,x_in,temp)
	  call rotate_z (lamdas,temp,x_out)
	else
	  call rotate_z (iverse*lamdas,x_in,temp)
	  call rotate_x (iverse*epsilon,temp,x_out)
	end if
c
	return
	end
