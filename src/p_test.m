function [x,optvalue,tpresolve,tphaseI,tphaseII,iter1,iter2,error_solution,error_outcome] = p_test(A,b,c)
S = size(A);
fprintf("Calculating the solution with standard MATLAB function...\n");
[~,fval] = linprog(c,-speye(S(1,2)),zeros(S(1,2),1),A,b);
fprintf("Finish solving the problem with standard MATLAB function!\n");
fprintf("Starting to solve the problem with my primal simplex method...\n");
[x,optvalue,tpresolve,tphaseI,tphaseII,iter1,iter2] = p_simplex(A,c,b);
error_solution = norm(A * x - b) / norm(x);
error_outcome = abs(optvalue - fval) / abs(fval);
end