c
c	transformation from GSE to GEO coordinates
c
c	itime(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
c	itime(2)=msec		miliseconds of day in UT
c	coordinates x_in and x_out are Cartesian and real*4
c
	subroutine gse_to_geo (itime,x_in,x_out)
c
	real*4 x_in(3),x_out(3),temp(3)
	integer*4 itime(2)
c
	call t2(itime,x_in,temp,-1)
	call t1(itime,temp,x_out,1)
c
	return
	end
