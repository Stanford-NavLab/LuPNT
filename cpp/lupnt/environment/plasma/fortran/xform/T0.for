c
c	function subroutine for computing time in Julian centuries
c
c	itime(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
c	itime(2)=msec		miliseconds of day in UT
c
	function t0(itime,iyr,iday,ut)
	integer*4 itime(2),iyr,iday
	real ut,t0,rmjd,fracday
c
c get time in Julian centuries
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
	T0=(RMJD-51544.5)/36525.0
c
	return
	end
