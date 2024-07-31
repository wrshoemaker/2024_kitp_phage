
import config
import numpy
import pylab as p
# Definition of parameters
a = 1.
b = 0.1
c = 1.5
d = 0.1
def dX_dt(X, t=0):
    """ Return the growth rate of fox and rabbit populations. """
    return numpy.array([ a*X[0] -   b*X[0]*X[1] ,
                  -c*X[1] + d*b*X[0]*X[1] ])



#!python
X_f0 = numpy.array([     0. ,  0.])
X_f1 = numpy.array([ c/(d*b), a/b])
all(dX_dt(X_f0) == numpy.zeros(2) ) and all(dX_dt(X_f1) == numpy.zeros(2)) # => True


#!python
def d2X_dt2(X, t=0):
    """ Return the Jacobian matrix evaluated in X. """
    return numpy.array([[a -b*X[1],   -b*X[0]     ],
                  [b*d*X[1] ,   -c +b*d*X[0]] ])


A_f0 = d2X_dt2(X_f0)


A_f1 = d2X_dt2(X_f1)                    # >>> array([[ 0.  , -2.  ],
                                        #            [ 0.75,  0.  ]])
# whose eigenvalues are +/- sqrt(c*a).j:
lambda1, lambda2 = numpy.linalg.eigvals(A_f1) # >>> (1.22474j, -1.22474j)
# They are imaginary numbers. The fox and rabbit populations are periodic as follows from further
# analysis. Their period is given by:
T_f1 = 2*numpy.pi/abs(lambda1)  


from scipy import integrate
t = numpy.linspace(0, 15,  1000)              # time
X0 = numpy.array([10, 5])                     # initials conditions: 10 rabbits and 5 foxes
X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
infodict['message']   



rabbits, foxes = X.T
f1 = p.figure()
p.plot(t, rabbits, 'r-', label='Rabbits')
p.plot(t, foxes  , 'b-', label='Foxes')
p.grid()
p.legend(loc='best')
p.xlabel('time')
p.ylabel('population')
p.title('Evolution of fox and rabbit populations')
f1.savefig('%srabbits_and_foxes_1.png' % config.analysis_directory)