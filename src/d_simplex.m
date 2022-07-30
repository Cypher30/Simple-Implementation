function [x,optvalue,iter,t] = d_simplex(A,b,c,I)
tic;
tol = 10^(-9);
A = full(A);
S = size(A);
B = A(:,I);
invB = inv(B);
x = zeros(S(1,2),1);
x(I) = B \ b;
cb = c(I);
iter = 0;
while 1
    p = cb' * invB;
    for i = 1:S(1,1)
        if x(I(i)) < -tol
            break;
        end
    end
    if x(I(i)) >= -tol
        fprintf("Optimal solution found!\n");
        break;
    end
    u = invB(i,:) * A;
    for j = 1:S(1,2)
        if u(1,j) < -tol
            break;
        end
    end
    if u(1,j) >= -tol
        fprintf("Primal infeasible!\n");
        x = [];
        optvalue = inf;
        return;
    else
        minJ = j;
        cbar = c(j) - p * A(:,j);
        beta = cbar / u(1,j);
        for k = j:S(1,2)
            if u(1,k) < -tol
                cbar = c(k) - p * A(:,k);
                temp = cbar / u(1,k);
                if (temp - beta > tol) || (cbar < 0)
                    if cbar < 0
                        beta = 0;
                        minJ = k;
                    else
                        minJ = k;
                        beta = temp;
                    end
                end
            end
        end
        v = invB * A(:,minJ);
        theta = x(I(i)) / v(i,1);
        x(minJ,1) = theta;
        x(I) = x(I) - theta * v;
        x(I(i)) = 0;
        I(i) = minJ;
        invB(i,:) = 1 / v(i,1) * invB(i,:);
        for k = 1:S(1,1)
            if k == i
                continue;
            else
                invB(k,:) = invB(k,:) - v(k,1) * invB(i,:);
            end
        end
        cb = c(I,1);
    end
    iter = iter + 1;
    if mod(iter,50) == 0
         invB = A(:,I) \ eye(S(1,1));
         x(I,1) = A(:,I) \ b;
    end
    if mod(iter,100) == 0
        fprintf("%d iterations, current cost: %f\n",iter,cb' * x(I,1));
    end
end
x(I) = A(:,I) \ b;
optvalue = c' * x;
t = toc;
end