Here are the interface of all my functions:

d_simplex
[x,optvalue,iter,t] = d_simplex(A,b,c,I)
input:
	A	the constraint matrix(full row rank)
	b	the constraint vector
	c	the object function
	I	the index of initial basis
output:
	x	the optimal solution
	optvalue	the optimal value
	iter	number of iterations
	t	time consumed in the process

d_test
[x,optvalue,iter,t,error_solution,error_outcome] = d_test(A,b,c,I)
input:
	A	the constraint matrix(full row rank)
	b	the constraint vector
	c	the object function
	I	the index of initial basis
output:
	x	the optimal solution
	optvalue	the optimal value
	iter	number of iterations
	t	time consumed in the process
	error_solution	||Ax - b|| / ||x||
	error_outcome	|optvalue - fval| / |fval| where optvalue and fval are the optimal value attained by my implementation and the standard function

p_simplex
[X,optvalue,tpresolve,tphaseI,tphaseII,iter1,iter2] = p_simplex(A,C,b)
input:
	A	the constraint matrix
	C	the object function
	b	the constraint vector
output:
	X	the optimal solution	
	optvalue	the optimal value
	tpresolve time consumed in the presolving
	tphaseI	time consumed in the phaseI
	tphaseII	time consumed in the phaseII
	iter1	number of iterations in phaseI
	iter2	number of iterations in phaseII

p_test
[x,optvalue,tpresolve,tphaseI,tphaseII,iter1,iter2,error_solution,error_outcome] = p_test(A,b,c)
input:
	A	the constraint matrix
	b	the constraint vector
	c	the object function
output:
	x	the optimal solution	
	optvalue	the optimal value
	tpresolve time consumed in the presolving
	tphaseI	time consumed in the phaseI
	tphaseII	time consumed in the phaseII
	iter1	number of iterations in phaseI
	iter2	number of iterations in phaseII
	error_solution	||Ax - b|| / ||x||
	error_outcome	|optvalue - fval| / |fval| where optvalue and fval are the optimal value attained by my implementation and the standard function

phaseI
[I,iter] = phaseI(A,b)
input:
	A	the constraint matrix
	b	the constraint vector
output:
	I	Index of initial basis
	iter	number of iterations in the process

presolve
[A1,b1,Ix,valx,flag] = presolve(A,b)
input:
	A	the constraint matrix
	b	the constraint vector
output:
	A1	the presolved matrix
	b1	the presolved vector
	Ix	index of solved decision variables
	valx	the value of solved decision variables
	flag	indicator (0 means normal, 1 means the problem is infeasible)

revised_simplex
[x,optvalue,iter] = revised_simplex(A,b,c,I0)
input:
	A	the constraint matrix
	b	the constraint vector
	c	the object function
	I0	the index for initial basis
output:
	x	the optimal solution
	optvalue	the optimal value
	iter	number of iteration

revised_simplex_phaseI
[I0,iter] = revised_simplex_phaseI(A,b,c,I0)
input:
	A	the constraint matrix
	b	the constraint vector
	c	the object function
	I0	the index of initial basis
output:
	I0	the index of basis
	iter	number of iteration

test_instance
[A,b,c,I] = test_instance(m,n,d)
input:
	m	number of rows
	n	number of columns
	d	density
output:
	A	the constraint matrix
	b	the constraint vector
	c	the object function
	I	the index of initial basis for dual simplex

Using function "p_test" to test the primal simplex
Using function "d_test" to test the dual simplex. Notice that the input matrix A must has full row rank.
