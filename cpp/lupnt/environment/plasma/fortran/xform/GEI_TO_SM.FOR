c
c	transformation from GEI to SM coordinates
c
c	time(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
c	time(2)=msec		miliseconds of day in UT
c	coordinates x_in and x_out are Cartesian and real*4
c
	subroutine gei_to_sm (itime,x_in,x_out)
c
	real*4 x_in(3),x_out(3),temp1(3),temp2(3)
	integer*4 itime(2)
c
	call t2(itime,x_in,temp1,1)
	call t3(itime,temp1,temp2,1)
	call t4(itime,temp2,x_out,1)
c
	return
	end
