function [x,optvalue,iter,I0] = revised_simplex(A,b,c,I0)
tol = 10^(-5);
S = size(A);
x = zeros(S(1,2),1);
B = A(:,I0);
invB = inv(B);
x(I0,1) = B \ b;
cb = c(I0,1);
iter = 0;
while 1
    p = cb' * invB;
    for j = 1:S(1,2)
        if ~ismember(j,I0)
            temp = c(j,1) - p * A(:,j);
            if temp < -10^(-8)
                break;
            end
        end
    end
    if temp >= -10^(-8)
        break;
    end
    u = invB * A(:,j);
    for i = 1:S(1,1)
        if u(i,1) > tol
            break;
        end
    end
    if u(i,1) <= tol
        fprintf("optimal cost is -inf\n");
        optvalue = -inf;
        x = [];
        return;
    else
        minI = I0(i,1);
        theta = x(I0(i,1),1) / u(i,1);
        mink = i;
        for k = i:S(1,1)
            if (u(k,1) > tol && x(I0(k,1),1) / u(k,1) - theta < -10^(-10)) || (u(k,1) > tol && abs(x(I0(k,1),1) / u(k,1) - theta) < 10^(-10) && I0(k,1) < minI) || (u(k,1) > tol && x(I0(k,1),1) < 0 && I0(k,1) < minI)
                if x(I0(k,1),1) < 0
                    theta = 0;
                    minI = I0(k,1);
                    mink = k;
                else
                    theta = x(I0(k,1),1) / u(k,1);
                    minI = I0(k,1);
                    mink = k;
                end
            end
        end
        x(j,1) = theta;
        x(I0,1) = x(I0,1) - theta * u;
        x(minI,1) = 0;
        I0(mink,1) = j;
        invB(mink,:) = 1 / u(mink,1) * invB(mink,:);
        for i = 1:S(1,1)
            if i == mink
                continue;
            else
                invB(i,:) = invB(i,:) - u(i,1) * invB(mink,:);
            end
        end
        cb = c(I0,1);
    end
    iter = iter + 1;
    if mod(iter,50) == 0
         invB = A(:,I0) \ eye(S(1,1));
         x(I0,1) = A(:,I0) \ b;
    end
    if mod(iter,100) == 0
        fprintf("%d iterations, current cost: %f\n",iter,cb' * x(I0,1));
    end
end
x(I0,1) = A(:,I0) \ b;
optvalue = cb' * x(I0,1);
end