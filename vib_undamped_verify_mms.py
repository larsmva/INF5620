#!/usr/bin/env python

	
"""
	Task a)	
	The problem is given as
	u'' + (w**2)u = f(t) 
	where u is displacment from equilibrium, f is the external force on u.
	
	The compact discretized notation is written as 
	[D_tD_t u + ( w**2)u =f]^n ,
	which is equivalent to
	[ u(n+1) - 2u(n) + u(n-1) ]/dt**2 +(w**2)u(n) = f
	u(n+1) = 2u(n) - (w*dt)**2)u(n) -u(n-1) +f*dt**2
	Rearranging the expression we obtain
	u(n+1) = (2 - (w*dt)**2)u(n) -u(n-1) +f*dt**2
	
	The derivative of u(0) can be used to define an expression for u(-1).
	D_2t u'(0) = [u(1) - u(-1) ]/2dt = V  -> u(-1)= u(1) -2*dt*V  

	Thus the equation for u(1) becomes
	
	u(1) = 1/2 * ( 2 + (w*dt)**2 )I + f*dt**2) + dt*V

	with u(0) =I 

	Intial conditions restrict u_e(x,t)= u_e(t) = c*t+d
			
			u_e(t) = I =d
			u'_e(t)= V = c   

	D_tD_t u_e(t) = 0 ==> u_e(t)*w**2 = f(t) = w**2*I*t + w**2*V
	
	Task b)
	[DtDt(u)](n) = [DtDt(c*t +d)](n) = c*[DtDt(t)](n) +   [DtDt(d)](n)
	
	[DtDt(d)](n) = [d(n+1) -2*d(n) + d(n-1)]/dt**2 = 0 since d(n) = constant
	[DtDt(t)](n) = [t(n+1) -2*t(n) +t(n-1)]/dt**2 = 0 , since t(n+1) - t(n) = dt and t(n-1) - t(n) =-dt

	R = [DtDt(u)](n) - u(t(n)) + f(t(n)) =  f(t(n) ) - u(t(n)) = 0 ,since u''(t) = 0
.		
 
		  
"""
import nose.tools as nt
import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import scitools.std as st


def simple_visualize(t,u   ): 
	st.plot(t,   u,   'r--o',
	legend=['numerical', 'exact'],
     	xlabel='t',
     	ylabel='u',
     	title='theta=%g, dt=%g' % (theta, dt),
     	savefig='%s_%g.png' % (theta2name[theta], dt),
     	show=True)


	

V, t, I, w, dt , a ,b = sym.symbols('V t I w dt a b')  # global symbols
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = DtDt(u,dt) + w**2*u(t) - f
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""

    R =  ( u(t+dt)  - u(t) )/dt  + 0.5*dt*w**2*u(t)  -V - 0.5*dt*f
    return sym.simplify(R.subs(t,0))

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.sym.as_finite_diff(u(t).diff(t,2), dt)
    """  	
    return (u(t+dt)- 2*u(t) + u(t-dt) )/dt**2

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u
    print "Initial conditions u(0)=%s, u'(0)=%s:" %(u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))
    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)




def cubic():
	"""
	Task E) : 
	Since the truncation error for the first step is of the third order. Then the numerical cubic solution and the exact cubic solution will be different.
	Thus the answer is NO.

	u(0) = I --> d = I 
	u'(0) = V --> c = V
	"""
	main(lambda t : a*t**3 + b*t**2 + t*V  + I) 

def quadratic():
	"""
	Task d)
	the quadratic func 

	u = b*t**2 + c*t +d 
	must fullfill the conditions 
	u(0) = I  --> d = I 
	u'(0) = V  --> c = V
	 
	"""
	main(lambda t : b*t**2 + t*V + I)
 
def linear():
    main(lambda t: V*t + I)

def solver_vib_undamped( source, u_0, dudt_0, Nt , T, omega= math.pi*2 ):
	"""
	Task f)
	======== ==================================================
	Argument                Explanation
	======== ==================================================
	source   function of t.
	T        Final time in time mesh.
	Nt       Number of time steps.
	u_0 	 intial condition of u, i.e  u(0) = u_0
	dudt_0   intial condition of u', i.e u'(0) = dudt_0
	w        vibration frequency  u'' + w**2*u = f(t)
	"""
	dt = float(T/Nt)
	t_n = [ dt*i for i in range(Nt)]
	if source==None or source==0:
		source = lambda t: 0

	f_n = map(source,t_n)
	
	u_n = np.array( [ 0.0 for i in t_n])
	OMEGA= (omega*dt)**2

	u_n[0]= u_0
	u_n[1]= 0.5*2*u_n[0] - 0.5*OMEGA*u_n[0]  + 0.5*f_n[0]*dt**2 + dt*dudt_0
	for i in range(1,Nt-1):
		u_n[i+1] = 2*u_n[0]- OMEGA*u_n[i] - u_n[i-1] + f_n[i]*dt**2

	return u_n,t_n

def test_solver_vib_undamped_quadratic_solution():
	"""
	Task g)
	Correct source term for u''+w**2*u = f(t) with u(t) = b*t**2 + V*t + I is
	f(t) = 2*b +w**2(  b*t**2 + V*t + I )

	selecting 
	b = 2 ,V = 4 , I = 8 , w = math.pi*2	
	
	"""	
	source = lambda t : 4*(1+(math.pi)**2*(  2*t**2 + 4*t + 8 ))
	u_num,t = solver_vib_undamped(source, 8 , 4 , 300 , 1 )
	def exact_quadratic_solution(t ,b, V, I):
		return b*t**2 + V*t + I

	u_e = np.array([exact_quadratic_solution(i,2,4,8) for i in t])
	
	e_max =np.abs(u_e -u_num).max()
	print e_max
	#assert (e_max < 1e-14)
	nt.assert_almost_equal(e_max, 0, places=14)	
	
if __name__ == '__main__':
    test_solver_vib_undamped_quadratic_solution()
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--order', type=int, default=4)
    parser.add_argument('--form', type=str, default="linear",help="Specify the form of u, i.e linear, quadratic and cubic")
    inval = parser.parse_args() 
    if inval.form=="linear":
	linear()
    elif inval.form == "quadratic":
	quadratic()
    elif inval.form =="cubic":
	cubic()
    else :
	parser.help()







