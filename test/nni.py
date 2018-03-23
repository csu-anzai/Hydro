import Ngl as ngl
import numpy as np
import matplotlib.pyplot as plt

# declare a list
#x = np.random.rand(6,1)
#y = np.random.rand(6,1)
x = [.45, .47, .55, .57, .61, .63]
y = [.45, .53, .46, .56, .57, .48]

# declare function
def f(x,y): return np.power(x,3.) + np.power(y,2.)

# point to interpolate
xn = [.5]
yn = [.5]

# interpolate
fa = ngl.natgrid(x,y,f(x,y),xn,yn)

# print true value
s = 'The true value is ' + repr(f(.5,.5)) + ', and the \
approximation is ' + repr(fa)

print(s)

# plot
plt.plot(x,y)
plt.show()
