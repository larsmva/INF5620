#!/usr/bin/env python
# -*- coding: utf-8 -*-
def visulasize_pendulum(x,y,save_id=""):
	
	plt.plot(x,y,'r.')


	plt.axis('equal')
    	plt.xlabel(r'x ')
    	plt.ylabel(r'y ')
	#plt.legend([''], loc='upper left')
	plt.title("Pendulum Motion in Cartesian coordinates")
	plt.savefig("position_pendulum_%s.png"%str(save_id))
	#plt.show(block=True)
	
	plt.clf()
 

def visulasize_pendulum_angle(t,theta,save_id=""):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(t,theta)
	plt.xlabel(r't $ \quad \left[time\right]$',fontsize=16)
	plt.ylabel(r'$\theta \quad \left[degrees\right] $',fontsize=16)
	plt.xlim((t[0],t[-1]))
	plt.title("Pendulum angle at time t ")
	plt.savefig("pendulum_angle_%s.png"%str(save_id))

	#plt.show(block=True)
	

	plt.clf()


def visulasize_compared_pendulums(elastic_theta, nonelastic_theta,time ,save_id=""):
	"""
	Plots elastic vs non-elastic pendulum
	"""
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(time,elastic_theta,'r-', time ,nonelastic_theta,'b-')
	plt.legend(["Elastic pendulum","Non-elastic pendulum"])
	plt.xlabel(r't $ \quad \left[time\right]$',fontsize=16)
	plt.ylabel(r'$\theta \quad \left[degrees\right] $',fontsize=16)
	plt.xlim((time[0],time[-1]))
	plt.title("Comparison plot between elastic and non-elastic pendulum angles")
	plt.savefig("Elastic_n_nonelastic_angles_%s.png"%str(save_id))
	#plt.show(block=True)
	plt.clf()

		
	

def simulate(
    beta=0.9,                 # dimensionless parameter
    Theta=30,                 # initial angle in degrees
    epsilon=0,                # initial stretch of wire
    num_periods=6,            # simulate for num_periods
    time_steps_per_period=60, # time step resolution
    plot=True,                # make plots or not
    ):
	
	"""
	Task 1)
	The scaled equation is given as

	x'' = -beta/(1-beta)(1 - beta/|L|)x	
	y'' = -beta/(1-beta)(1 - beta/|L|)(y-1) -beta
	
	|L| = ( x**2 + ( y-1)**2)**0.5
	with inital conditions
	x(0) = (1 +epsilon)sin(Theta) 	, x'(0)=0
	y(0) = 1- (1 +epsilon)cos(Theta), y'(0)=0

	The length of the wire, denoted as L, is depenedent on timestep n.
	|L|= L --> L(t) --> L(x(t(n),y(t(n))) = L(n)

	The  discretized expression of second derivatives is
	u'' = [u(n+1) -2*u(n) + u(n-1) ]/dt**2
	
	This leads 
	[x(n+1) -2*x(n) + x(n-1) ]/dt**2 = -beta/(1-beta)(1-beta/L(n)x(n) 
	[y(n+1) -2*y(n) + y(n-1) ]/dt**2 = -beta/(1-beta)(1-beta/L(n))(y(n)-1) - beta 

	Rearranging these expressions yields
	x(n+1) =  [dt**2*beta/(1-beta)]*(1-beta/L(n))*x(n) +2*x(n) - x(n-1)
	y(n+1) =  [dt**2*beta/(1-beta)]*(1-beta/L(n))*(y(n)-1) - beta*dt**2 +2*y(n) - y(n-1)

	We use of the first derivative of x and y to find the values of x(-1) and y(-1). The first 	   derivatives are given as 
	D_2t = [u(n+1) - u(n-1)]/2*dt = 0  
	Thus 
		x(-1) = x(1)  and y(-1) ) y(1)

	Hence the equtions for x(1) and y(1)  becomes
	x(1) = -0.5*(dt**2)*beta/(1-beta)(1-beta/L(0))x(0) +x(0)
	y(1) = -0.5*(dt**2)*beta/(1-beta)(1-beta/L(0))(y(0)-1) +y(0) - 0.5*beta*dt**2
	

	"""

	P= np.pi*2
	time_steps = int(P*time_steps_per_period)
	time = num_periods*P
	dt = float(time / time_steps)
	t = [ i*dt for i in range(time_steps)] 



	
	x = [ 0.0 for i in t]
	y = [ 0.0 for i in t]
	x[0] = (1. +epsilon)*math.sin(math.radians(Theta))
	y[0] = 1.- (1. +epsilon)*math.cos(math.radians(Theta))

	""" 
	In the case of non-elastic pendulum , i.e beta --> 1, we ned to use L'Hoptials rule.
	lim ( beta --> 1) [1-beta/|L|]/[1-beta] = 1/|L|
	Thus 
	lim ( beta --> 1 )  [beta/(1-beta)(1-beta/L(n))] = beta/|L|
	
	With Theta=0 we have y-1= 0 . 

	So instead of handling both cases, it is sufficent to add a small constant to the value.

	"""

	dt2=dt**2
	
	if beta ==1 :
	  	beta -= 0.001

	if Theta == 0:
		Theta-= 0.0001

	
 	aux = beta/(1-beta)
	
	L = lambda x,y : (x**2+(y-1.0)**2)**0.5	
	
	dt2=dt**2
	Linv = beta/L(x[0],y[0])

	x[1] = -0.5*dt2*aux*(1- Linv)*x[0]      +x[0]
	y[1] = -0.5*dt2*aux*(1- Linv)*(y[0]-1)  +y[0] - 0.5*beta*dt2

	for n in range(1,time_steps-1) : 
		Linv = beta/L(x[n],y[n])
		x[n+1] = -dt2*aux*(1.- Linv)*x[n] + 2.0*x[n] - x[n-1]
		y[n+1] = -dt2*aux*(1.- Linv)*(y[n]-1)  + 2.0*y[n] - y[n-1] - beta*dt2

	
	theta_func = lambda x,y : math.degrees(math.atan(x/(1-y)))  
	theta = map(theta_func,x,y)
	

	if plot : 
		visulasize_pendulum(x,y,save_id="beta%f_theta%f_epsilon%f"%(beta,Theta,epsilon))
		visulasize_pendulum_angle(t,theta,save_id="beta%.3f_theta%.3f_epsilon%.3f"%(beta,Theta,epsilon))
		if ( Theta < 10 ) :
			x_non,y_non,theta_non=simulate(1,Theta,epsilon, num_periods,  time_steps_per_period, plot=False)
			visulasize_compared_pendulums(theta,theta_non,t,save_id="beta%.3f_theta%.3f_epsilon%.3f"%(beta,Theta,epsilon))
		
			
	return x,y,theta


def test_stationary():
		beta=0.5
		x,y,theta = simulate(beta,0,plot=False)
	
		tol = 1e-13		
		assert (np.abs(x).max()  < tol  and np.abs(y).max() < tol and np.abs(theta).max()< tol)
		
		

def test_vertical_oscillations():
		"""  
			with Theta = 0 the initial conditions becomes
			x(0) = 0		, x'(0) = 0 
			y(0) = epsilon		, y'(0) = 0
			
			Given that 
			x'' = -beta/(1-beta)(1 - beta/|L|)x ,
			leads to x'' = 0. This combined with x'(0) = 0 and x(0) = 0.
			Thus x=0 for all x.  

			Therefore we have only the vertical component  
			y'' = -beta/(1-beta)(1 - beta/|L|)(y-1) -beta

			with |L| = ((y-1)**2)**0.5 = |y-1|
			y'' = 	-beta/(1-beta)(y-1) + beta**2/(1-beta)*( (y-1)/|y-1|)
			
			y'' = -y * beta/(1-beta) + beta/(1-beta) + beta**2/(1-beta)*( (y-1)/|y-1|) - beta

			Given that y < 1 , since it is scaled for (0,1) , thus  y ~ 0 
			This gives (y-1)/|y-1| = -1 , 

			y'' + y * beta/(1-beta) = + beta/(1-beta)  -beta**2/(1-beta) - beta
	
			y'' + y * beta/(1-beta) = beta/(1-beta) ( 1 - beta - (1- beta) ) = 0
			
			The differential vibration equation is given as

			y'' + y * w**2
						
			Thus with w**2 = beta/(1-beta) =>  w = [ beta/(1-beta)]**0.5
							 
			we obtain vertical oscillations / vibrations. The general solution for the vibration eqaution is 

			y = A cos(w*t) +B*sin(w*t)
			with the inital consitions 
			y(0) = - epsilon --> A = -epsilon
			y'(0) = 0 --> B =0
			We obtain the exact solution
			y = - epsilon*cos(w*t)
		"""
		time_steps_per_period =60000
		num_periods=2
		beta,epsilon = 0.8,0.1
		x,y,theta = simulate(beta,0,epsilon,num_periods,time_steps_per_period,plot=False )
	
		P= np.pi*2
		time_steps =int( P*time_steps_per_period)
		time = num_periods*P
		dt = float(time / time_steps)

		t = [ i*dt for i in range(time_steps)] 

		omega = (beta/(1-beta))**0.5

		def vertical_vib_exact_solution(t,epsilon,omega):
			return -epsilon*math.cos(omega*t)

		ye = np.array([ vertical_vib_exact_solution(i,epsilon,omega) for i in t])
		error_max  = np.abs(np.array(y) -ye).max()
		assert (error_max < 1e-7)	
	


def demo ( beta,Theta):
		"""
		Function demo takes 2 arguments 

		---------------------------------

		beta :

		Theta: The intial angle of the pendelum. 

		"""
	
		x,y,theta = simulate(beta,Theta,num_periods=3,time_steps_per_period=600)
 

		




if __name__ =='__main__':
    import argparse
    import numpy as np
    import matplotlib.pyplot as plt
    import math

    parser = argparse.ArgumentParser()
    parser.add_argument('--beta', type=float, default=0.9, help = " Float value of the scaled variabel beta.")
    parser.add_argument('--Theta', type=float, default=30.0, help = " Float value the inital angle displacment in degrees")
    parser.add_argument('--epsilon', type=float, default=0.0, help = "Float value of the radial displacement of the pendulum")
    parser.add_argument('--num_periods', type=int, default=6, help = "Integer of number of periods")
    parser.add_argument('--time_steps_per_period', type=int, default=60, help = "Integer of number of time steps pr. period")
    parser.add_argument('--plot', type=bool, default=True, help="Boolean value, True will display plots, while False will not" )
    parser.add_argument('--operation', type=str, default="simulate", help="  Options : demo , vertical_vib, stationary,simulate")
    inVal = parser.parse_args() 
    if inVal.operation=="demo":
	demo(inVal.beta,inVal.Theta,)
    elif inVal.operation=="vertical_vib":
	test_vertical_oscillations()
    elif inVal.operation=="stationary":
	test_stationary()
    elif inVal.operation=="simulate":
	simulate(inVal.beta,inVal.Theta,inVal.epsilon,inVal.num_periods, inVal.time_steps_per_period,inVal.plot)
    else :
	parser.help()
    









