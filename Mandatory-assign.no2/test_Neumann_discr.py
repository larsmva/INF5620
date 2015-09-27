#!/usr/bin/env python
from Neumann_discr import *

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
		import nose.tools as nt
		e = abs(u-exact_solution(x,t[n])).max()
		nt.assert_almost_equal(e, 0, places=13)

	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=sum_approximation)
	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=centered_difference)
	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=one_sided_difference)
	solve_wave_eqn_with_variable_velocity(q,f,I, L,T,Nx,u_exact=exact_solution ,user_action=assert_if_error,neumann=shifted_domain)

def test_convergence_rates_():
	
	import math
	import numpy as np
	import sympy
	from sympy import cos
	from sympy import sin
	from sympy.utilities.lambdify import lambdify
	import nose.tools as nt

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
		import nose.tools as nt
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

	

