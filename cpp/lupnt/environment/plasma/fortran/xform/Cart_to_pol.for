c
	subroutine cart_to_pol(x_in,lat,long,radial)
	parameter (pi=3.1415927)
	real lat,long,radial,x_in(3),row
c
	row=sqrt(x_in(1)*x_in(1)+x_in(2)*x_in(2))
	radial=sqrt(row*row+x_in(3)*x_in(3))
	lat=atan2(x_in(3),row)
	long=atan2(x_in(2),x_in(1))
	long=long + (1.0-sign(1.0,long))*pi
c
	return
	end
