c
c	transformation matrix t5=<phi-90,Y>*<lamda,Z>
c
c	iverse = 1 forward transform
c	       =-1 inverse transform
c
	subroutine t5 (itime,x_in,x_out,iverse)
	real x_in(3),x_out(3),phi,lamda,rmjd,ut
	real temp(3),fracday
	integer*4 itime(2),iverse,iyr,iday
	parameter (pid2=3.1415927/2.0,degrad=3.1415927/180.0)
c
	IYR=itime(1)/1000
	IDAY=itime(1)-IYR*1000
	ut=itime(2)/3600000.0
	fracday=ut/24.0
c
c	RMJD=44238.0+FLOAT(IYR-1980)*365.+
c     &		FLOAT((IYR-1977)/4)+FLOAT(IDAY)
c  rmjd = modified Julian day = time measured in days from 00:00UT November 17, 1858
	rmjd=45.0+float(iyr-1859)*365.0 + (float((iyr-1861)/4)+1.0)
     &		+ float(iday) -1.0 + fracday
c
	factor=(rmjd-46066.0)/365.25
	phi=(78.8+4.283e-2*factor)*degrad
	lamda=(289.1-1.413e-2*factor)*degrad
c
	if (iverse.eq.1) then
	  call rotate_z(lamda,x_in,temp)
	  call rotate_y(phi-pid2,temp,x_out)
	else
	  call rotate_y(iverse*(phi-pid2),x_in,temp)
	  call rotate_z(iverse*lamda,temp,x_out)
	end if
c
	return
	end
