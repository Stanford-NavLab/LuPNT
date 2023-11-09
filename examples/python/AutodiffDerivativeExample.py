import pylupnt as pnt
import sys

sys.path.append("pylupntutil")
import pntautodiff as ad

def myfunc(x):
    y =x + x*x + 1/x + pnt.log(x)   # pnt dual
    return y

def myfunc2(x, y, z):
    return 1 + x + y + z + x*y + y*z + x*z + x*y*z + pnt.exp(x/y + y/z)


print("derivatives of single value function")
print("--------------------------------------")
x = pnt.real(2.0)
y = myfunc(x)
print("y = ", y)
dydx = pnt.derivative(myfunc, x)
print("dydx = ", dydx)
print(" ")

print("derivatives of multi-value function")
print("--------------------------------------")
x = pnt.real(1.0)
y = pnt.real(2.0)
z = pnt.real(3.0)
z = myfunc2(x, y, z)
print("z = ", z)
dzdx = pnt.derivative3(myfunc2, x, y, z, wrt=0)  # 3 indicates the number of arguments
dzdy = pnt.derivative3(myfunc2, x, y, z, wrt=1)
dzdz = pnt.derivative3(myfunc2, x, y, z, wrt=2)
print("dz/dx = ", dzdx)
print("dz/dy = ", dzdy)
print("dz/dz = ", dzdz)

