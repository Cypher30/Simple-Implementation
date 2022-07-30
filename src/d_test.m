function [x,optvalue,iter,t,error_solution,error_outcome] = d_test(A,b,c,I)
S = size(A);
fprintf("Calculating the solution with standard MATLAB function...\n");
[~,fval] = linprog(c,-speye(S(1,2)),zeros(S(1,2),1),A,b);
fprintf("Finish solving the problem with standard MATLAB function!\n");
fprintf("Starting to solve the problem with my dual simplex method...\n");
[x,optvalue,iter,t] = d_simplex(A,b,c,I);
error_solution = norm(A * x - b) / norm(x);
error_outcome = abs(optvalue - fval) / abs(fval);
end