#!/usr/bin/env python
"""

The wave equation is given as

	u_tt = c**2 * u_xx +f 

with c as the wave velocity.

For dynamic wave velocity we have that

	u_tt = dx ( q(x)*u_x) + f(t,x)


With center discretization 

	[DtDt u](n,i) = [ Dx (q Dx u )](n,i) +f(n,i)
	
	we obtain

	u(n+1,i)-2*u(n,i)+u(n-1,i) = C**2 ( q(n,i+1/2)(u(n,i+1) -u(n,i) ) - q(n,i-1/2) ( u(n,i) - u(n,i-1) )   *1*

	with C = (dt/dx).


The estiametion of q(n,i+1/2) and q(n,i-1/2) can be done with the mean and here is 3 ways to calculate the mean:
	
	arithmetic mean q(n,i+1/2) = 0.5 *( q(n,i+1) + q(n,i)) = bar(q)

	harmonic mean q(n,i+1/2) = 2*( 1/q(n,i+1) + 1/q(n,i))**-1

	geometric mean q(n,i+1/2) = (q(n,i+1)*q(n,i))**0.5


We will use the arithmetic mean, thus the explicit scheme becomes

	u(n+1,i) = 2*u(n,i)- u(n-1,i) + C**2 ( 0.5 *( q(n,i+1) + q(n,i)) (u(n,i+1) -u(n,i) ) -  0.5 *( q(n,i) + q(n,i-1)) ( u(n,i) - u(n,i-1) ) +f(n,i)

"""

"""
In the next function, we will have the same number of arguements.
These argumetns are explained here:

======== ==================================================
Argument                Explanation
======== ==================================================
u_1 	 The values of u for the previous time step.
q	 The squared wave velocity function
i 	 Spacial index
x	 Values in the spacial domain.
n	 Temporal index
eval_q	 Method for evaluation of q(i+1/2) and q(i-1/2)


Additional comments : 

"""
def centered_difference(u_1,q,i,x,n,eval_q):
		v = (1,-1)[int(i>0)]
		if n==1:
			return eval_q(q,x,i,v)*(u_1[i+v]-u_1[i]) 
		else : 
			return 2.0*eval_q(q,x,i,v)*(u_1[i+v]-u_1[i])

def sum_approximation(u_1,q,i,x,n,eval_q):
		v = (1,-1)[int(i>0)]
		if n==1:
			return eval_q(q,x,i,0)*(u_1[i+v]-u_1[i]) 
		else : 
			return 2*eval_q(q,x,i,0)*(u_1[i+v]-u_1[i]) 

def one_sided_difference(u_1,q,i,x,n,eval_q):
		v = (1,-1)[int(i>0)]
		if n==1:
			return 0.5*eval_q(q,x,i,0)*(u_1[i+v]-u_1[i]) 
		else : 
			return eval_q(q,x,i,0)*(u_1[i+v]-u_1[i]) 
			
def shifted_domain(u_1,q,i,x,n,eval_q):
		v = (1,-1)[int(i>0)]
		if n==1:
			return  0.5*eval_q(q,x,i,v)*(u_1[i+v]-u_1[i]) 
		else : 
			return  eval_q(q,x,i,v)*(u_1[i+v]-u_1[i]) 



def test_constant():
	import math
	import numpy as np


	exact_solution = lambda x,t : 2.345
	q = lambda x : 4.245
	f =  lambda x,t : 0
	I = lambda x : 2.345
	L = 2
	T=1
	Nx = 50  

	def assert_if_error(u,x,t,n):
		e = abs(u-exact_solution(x,t[n])).max()
		nt.assert_almost_equal(e, 0, places=13)

	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=sum_approximation)
	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=centered_difference)
	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=one_sided_difference)
	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=shifted_domain)


def error_of_methods_with_symmteric_boundary():
	"""
	Complie test with returning error
	"""

	import math
	import numpy as np
	import sympy
	from sympy import cos
	from sympy import sin
	from sympy.utilities.lambdify import lambdify
	

	w,x, L ,t= sympy.symbols("w x L t")
	pi = math.pi	

	u = lambda x,t : cos(w*t)*cos(pi*x/L)
	q = lambda x : 1 + cos(pi*x/L)

	def source_term(u,q):
		return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

	w=1
	L = 8
	T=4
	Nx=40

	f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
	q = sympy.lambdify(x,q(x) ,'numpy')
	I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
	u = sympy.lambdify((x,t), u(x,t),'numpy')

	
	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=sum_approximation)[1]
	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=centered_difference)[1]
	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=centered_difference)[1]
	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=shifted_domain)[1]

	

def error_of_methods():
	"""
	Complie test with returning error
	"""
	import math
	import numpy as np
	import sympy
	from sympy import cos
	from sympy import sin
	from sympy.utilities.lambdify import lambdify
	
	w,x, L ,t= sympy.symbols("w x L t")
	pi = math.pi	

	u = lambda x,t : cos(w*t)*cos(pi*x/L)
	q = lambda x : 1 + ( x -L/2.0)**4

	def source_term(u,q):
		return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

	w=1
	L = 2
	T=10
	Nx=40

	f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
	q = sympy.lambdify(x,q(x) ,'numpy')
	I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
	u = sympy.lambdify((x,t), u(x,t),'numpy')

	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=sum_approximation)[1]
	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=centered_difference)[1]
	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=centered_difference)[1]
	print solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=u,neumann=shifted_domain)[1]



def solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=None,user_action=None,neumann=None):
	"""
	======== ==================================================
	Argument                Explanation
	======== ==================================================
	q        Squared wave velocity, either function or array of
		 values at mesh points
	f 	 the source term
	I 	 the inital value at t=0
	L 	 length of the spacial domain [0,L]
	T 	 Final time of the temporal domain [0,T]
	Nx 	 Number of discretization points in the spacaial 
		 domain. 
	u_exact  The exact solution.
    user_action  Function that is called at certain points in the 
		 solver.
       neumann   The specific Neumann boundary conditions in the  
		 spacial domain.
	"""
	import numpy as np
	import inspect
	x ,dx= np.linspace(0, L, Nx, retstep=True)
	
	
	def stability(q,x,beta):
		return beta/max(map(q,x))**0.5

	dt = dx * stability(q,x,0.9)

	t ,dt= np.linspace(0, T,int(T/dt) ,retstep=True)

	#q = np.array(map(q,x))
	if f==None:
		f = lambda x,t:0


	"""
	The following if-sentence is to establish if q is a function or 
	a collection of values at mesh points. Then selecting the accurate 
	method to evauate q at q(i+1/2) and q(i-1/2)

	"""
	if inspect.isfunction(q):
		
		def half_step(q,x,i,offset=0):
			return q(x[i]+offset*dx/2)	
		eval_q=half_step
	else:
		if len(q) == len(x):
			def mean(q,x,i,offset=0):
				return 0.5*(q[i] +q[i+offset])
			eval_q = mean 
		else:
			return
 
	

	u   = np.zeros(Nx)
	if inspect.isfunction(I):
		u_1 = np.array(map(I,x))
	else:
		u_1 = I
	u_2 = np.zeros(Nx)

	C2= (dt/dx)**2
	dt2= dt**2
	#q = np.array(map(q,x))
	
	Inp = np.array(range(1,Nx-1)) # Internal points 
	
	u[0]  = u_1[0] + C2*neumann(u_1,q,   0,x,1,eval_q)     + 0.5*dt2*f(x[0],t[0])
	u[-1] = u_1[-1]+ C2*neumann(u_1,q,Nx-1,x,1,eval_q) + 0.5*dt2*f(x[-1],t[0])

	u[Inp] = u_1[Inp] +\
			0.5*C2*(eval_q(q,x,Inp,1)*(u_1[Inp+1] - u_1[Inp]) - \
				eval_q(q,x,Inp,-1)*(u_1[Inp] - u_1[Inp-1])) + \
 								   0.5*dt2*f(x[Inp], t[0])

	if user_action!=None:
		user_action(u,x,t,0)

	u_2, u_1, u = u_1, u, u_2

	
	for n in range(1,len(t)):
		u[Inp] = - u_2[Inp] + 2*u_1[Inp] + \
				+C2*( eval_q(q,x,Inp,1)*(u_1[Inp+1] - u_1[Inp]) - \
					eval_q(q,x,Inp,-1)*(u_1[Inp] - u_1[Inp-1]) ) + \
									   dt2*f(x[Inp],t[n])

		u[-1] = 2.0*u_1[-1]- u_2[-1] + C2*neumann(u_1,q,Nx-1,x,n,eval_q)  + dt2*f(x[-1],t[n])
		u[0]  = 2.0*u_1[0] - u_2[0]  + C2*neumann(u_1,q, 0  ,x,n,eval_q)  + dt2*f(x[0],t[n])

		if user_action!=None:
			user_action(u,x,t,n)
		u_2, u_1, u = u_1, u, u_2
		
		
	e = abs(u_1-u_exact(x,t[-1]))	
	e  =sum(e**2)*dx

	return u ,e, dt


def test_convergence_rates_symmteric_boundary():

	import math
	import numpy as np
	import sympy
	from sympy import cos
	from sympy import sin
	from sympy.utilities.lambdify import lambdify
	

	w,x, L ,t= sympy.symbols("w x L t")
	pi = math.pi	

	u = lambda x,t : cos(w*t)*cos(pi*x/L)
	q = lambda x : 1 + cos(pi*x/L)

	def source_term(u,q):
		return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

	w=1
	L = 2
	T=2
	Nx=50

	f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
	q = sympy.lambdify(x,q(x) ,'numpy')
	I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
	exact_solution = sympy.lambdify((x,t), u(x,t),'numpy')

	r_exact=2.0	
	def assert_convergence(r):
		tol= 2.5
		for i in range(len(r)):
			diff = abs(r_exact -r[i])
			assert(  diff <= (tol+0.5) ) 
			tol=diff



	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=sum_approximation,user_action=assert_convergence)
	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=centered_difference,user_action=assert_convergence)
	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=one_sided_difference,user_action=assert_convergence)
	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=shifted_domain,user_action=assert_convergence)


def test_convergence_rates_():
	import math
	import numpy as np
	import sympy
	from sympy import cos
	from sympy import sin
	from sympy.utilities.lambdify import lambdify
	

	w,x, L ,t= sympy.symbols("w x L t")
	pi = math.pi	

	u = lambda x,t : cos(w*t)*cos(pi*x/L)
	q = lambda x : 1 + ( x -L/2.0)**4

	def source_term(u,q):
		return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

	w=1
	L=2
	T=2
	Nx=50

	r_exact=2.0	
	def assert_convergence(r):
		tol= 2.5
		for i in range(len(r)):
			diff = abs(r_exact -r[i])
			assert(  diff <= (tol+0.5)  )  # 0.5 because of oscillations 
			tol=diff
			


	f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
	q = sympy.lambdify(x,q(x) ,'numpy')
	I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
	exact_solution = sympy.lambdify((x,t), u(x,t),'numpy')


	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=sum_approximation,user_action=assert_convergence)
	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=centered_difference,user_action=assert_convergence)
	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=one_sided_difference,user_action=assert_convergence)
	compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=shifted_domain,user_action=assert_convergence)
	


def compute_convergence_rates(q,f,I, L,T,Nx,u,neumann=None,user_action=None):
	"""
	======== ==================================================
	Argument                Explanation
	======== ==================================================
	q        Squared wave velocity, either function or array of
		 values at mesh points
	f 	 the source term
	I 	 the inital value at t=0
	L 	 length of the spacial domain [0,L]
	T 	 Final time of the temporal domain [0,T]
	Nx 	 Number of discretization points in the spacaial 
		 domain. 
	u        The exact solution, not optional
    user_action  Function that is called at certain points in the 
		 solver.
       neumann   The specific Neumann boundary conditions in the  
		 spacial domain.
	"""
	if u==None: raise ValueError
	E = []
	h=[]
	for nx in [ Nx*i for i in range(1,10)]:
	
		u_num, e , dt =solve_wave_eqn_with_variable_velocity(q,f,I, L,T,nx,u_exact=u,neumann=neumann)
		E.append(e)
		h.append(dt)
	
	
	r=  convegence_rates(h,E)

	if user_action!=None:
		user_action(r)

	return r

def convegence_rates(h,E):
	from math import log
	return [log(E[i]/E[i-1])/log(h[i]/h[i-1]) for i in range(1, len(h))]

if __name__ == '__main__':
    import argparse 
    import math
    import numpy as np
    import sympy
    from sympy import cos
    from sympy import sin
    from sympy.utilities.lambdify import lambdify
    parser = argparse.ArgumentParser()
  
    parser.add_argument('--T', type=float, default=2.0 ,help="Lenght of the temporal domain ") 
    parser.add_argument('--L', type=float, default=2.0 ,help="Lenght of the spacial  domain ")
    parser.add_argument('--Nx', type=int, default=50 ,help="Number of mesh points in the spacial domain ")
    parser.add_argument('--task', type=str, default="13A" ,help="Options : all, 13A, 13B, 13C, 13D, test-scalar, error-methods,test-convergence ")
    parser.add_argument('--dx', type=float, default=None ,help="The discretaztion length in the spacial domain, superseeds the input of Nx ")
 
    Z = parser.parse_args() 
    if Z.dx!=None:
	Z.Nx = int(round(N.L/N.dx))

    task_all=False
    if Z.task=="all":
	task_all=True
	
    if Z.task=="13A" or task_all:
	def  task_13A(in_L,in_T,Nx):

	

		w,x, L ,t= sympy.symbols("w x L t")
		pi = math.pi	

		u = lambda x,t : cos(w*t)*cos(pi*x/L)
		q = lambda x : 1 + ( x -L/2.0)**4

		def source_term(u,q):
			return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

		w=1
		L,T = in_L,in_T
		f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
		q = sympy.lambdify(x,q(x) ,'numpy')
		I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
		exact_solution = sympy.lambdify((x,t), u(x,t),'numpy')

		r = compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=sum_approximation)
		print "====> Convergences rate for approx q(i+1/2) + q(i-1/2)=2q(i) :" , r[-1] 
		r=compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=centered_difference)
		print "====> Convergences rate for approx  q(i+1/2)= q(i-1/2) :" , r[-1] 
	task_13A(Z.L,Z.T,Z.Nx)

    if Z.task=="13B" or task_all:
    	def task_13B(in_L,in_T,Nx):

		w,x, L ,t= sympy.symbols("w x L t")
		pi = math.pi	

		u = lambda x,t : cos(w*t)*cos(pi*x/L)
		q = lambda x : 1 + cos(pi*x/L)

		def source_term(u,q):
			return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

		w=1
		L,T = in_L,in_T

		f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
		q = sympy.lambdify(x,q(x) ,'numpy')
		I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
		exact_solution = sympy.lambdify((x,t), u(x,t),'numpy')

		r = compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=centered_difference)
		print "====> Convergences rate for approx  q(i+1/2)= q(i-1/2) :" , r[-1] 
	task_13B(Z.L,Z.T,Z.Nx)
    if Z.task=="13C" or task_all:
	def task_13C(in_L,in_T,Nx):

		w,x, L ,t= sympy.symbols("w x L t")
		pi = math.pi	

		u = lambda x,t : cos(w*t)*cos(pi*x/L)
		q = lambda x : 1 + cos(pi*x/L)

		def source_term(u,q):
			return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

		w=1
		L,T = in_L,in_T

		f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
		q = sympy.lambdify(x,q(x) ,'numpy')
		I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
		exact_solution = sympy.lambdify((x,t), u(x,t),'numpy')

		r = compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=one_sided_difference)
		print "====> Convergences rate for one sided difference :" , r[-1] 
	task_13C(Z.L,Z.T,Z.Nx)
    if Z.task=="13D" or task_all:
	def task_13D(in_L,in_T,Nx):

		w,x, L ,t= sympy.symbols("w x L t")
		pi = math.pi	

		u = lambda x,t : cos(w*t)*cos(pi*x/L)
		q = lambda x : 1 + cos(pi*x/L)

		def source_term(u,q):
			return sympy.simplify(u(x,t).diff(t, t) - (u(x,t).diff(x)*q(x)).diff(x) )

		w=1
		L,T = in_L,in_T

		f =  sympy.lambdify((x,t), source_term(u,q)  ,'numpy')
		q = sympy.lambdify(x,q(x) ,'numpy')
		I = sympy.lambdify(x, u(x,t).subs(t,0),'numpy')
		exact_solution = sympy.lambdify((x,t), u(x,t),'numpy')

		r = compute_convergence_rates(q,f,I, L,T,Nx,exact_solution,neumann=shifted_domain)
		print "====> Convergences rate for domain shift:" , r[-1] 
	
	task_13D(Z.L,Z.T,Z.Nx)

    if Z.task=="test-scalar":
	test_constant()
    if Z.task=="error-methods":
	error_of_methods()
	error_of_methods_with_symmteric_boundary()
    if Z.task=="test-convergence":
    	test_convergence_rates_()
	test_convergence_rates_symmteric_boundary()
    
	
