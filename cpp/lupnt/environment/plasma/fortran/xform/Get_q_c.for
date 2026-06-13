c
c	subroutine for computing the dipolar axis in GSE coordinates
c
	subroutine get_q_c (itime,q_c)
	real q_g(3),q_c(3),ut,rmjd,factor,phi,lamda,temp(3)
	integer*4 itime(2)
	parameter (degrad=3.1415927/180.0)
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
	call pol_to_cart(phi,lamda,1.0,q_g)
	call t1(itime,q_g,temp,-1)
	call t2(itime,temp,q_c,1)
c
	return
	end
