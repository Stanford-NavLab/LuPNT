
c
c	transformation from GEI to GEO coordinates
c
c	time(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
c	time(2)=msec		miliseconds of day in UT
c	coordinates x_in and x_out are Cartesian and real*4
c
	subroutine gei_to_geo (itime,x_in,x_out)
c
	real*4 x_in(3),x_out(3)
	integer*4 itime(2)
c
	call t1(itime,x_in,x_out,1)
c
	return
	end
