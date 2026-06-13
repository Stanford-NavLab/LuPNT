c
c	Function varies from 0 to 1 to within 0.1%
c	as the arguement x passes from (a-da) to (a+da).
c	So "da" is the range over which the function turns on.
c
	function switchon(x,a,da)
c
	real x,a,da,c,switchon
c
	c=3.4534/da
	switchon=tanh(c*(x-a))/2+0.5
	return
	end
