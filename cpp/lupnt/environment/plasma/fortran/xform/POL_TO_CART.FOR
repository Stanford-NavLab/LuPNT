c
	subroutine pol_to_cart(lat,long,radial,x_out)
	real lat,long,x_out(3)
c
	coslat=cos(lat)
c
	x_out(1)=radial*coslat*cos(long)
	x_out(2)=radial*coslat*sin(long)
	x_out(3)=radial*sin(lat)
c
	return
	end
