from sympy import *

# 10.6 and 10.7 solved here

rho = var("rho")
x,y,z = var("x,y,z")
a, h = var("a, h")
m = var("m")

task = "10.6"
task = "10.7"

if task == "10.6":
    i = y*z
elif task == "10.7":
    i = x*y

i1 = integrate( i, (y, 0, a-x) )
i2 = integrate( i1, (x, 0, a) )
i3 = integrate( i2, (z, 0, h) )
I = rho*i3

if task == "10.6":
    pprint("\n10.6: Iyz:")
    pprint(I)

elif task == "10.7":
    pprint("\n10.7: Ixy:")
    pprint(I)
    I = I.subs(rho/2*a*a*h, m)
    pprint("\n")
    pprint(I)
