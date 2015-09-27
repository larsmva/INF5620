from Neumann_discr import *

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


    if Z.task=="error-methods":
	error_of_methods()
	error_of_methods_with_symmteric_boundary()
    if Z.task=="test-convergence":
	import test_Neumann_discr as tst

    	tst.test_convergence_rates_()
	tst.test_convergence_rates_symmteric_boundary()

    if Z.task=="test-scalar":
	import test_Neumann_discr as tst
	tst.test_constant()
