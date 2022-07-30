function [I,iter] = phaseI(A,b)
tol = 10^(-10);
S = size(A);
for i = 1:S(1,1)
    if b(i,1) < 0
        A(i,:) = -A(i,:);
        b(i,:) = -b(i,:);
    end
end
A1 = zeros(S(1,1),S(1,1) + S(1,2));
A1(:,1:S(1,2)) = A;
A1(:,S(1,2) + 1:S(1,2) + S(1,1)) = eye(S(1,1));
I = S(1,2) + 1:S(1,2) + S(1,1);
I = I';
c = zeros(S(1,1) + S(1,2),1);
c(S(1,2) + 1:S(1,1) + S(1,2),1) = 1;
[I,iter] = revised_simplex_phaseI(A1,b,c,I);
B = A1(:,I);
invB = B \ eye(S(1,1));
for k = 1:length(I)
    if I(k,1) > S(1,2)
        for j = 1:S(1,2)
            if abs(invB(k,:) * A(:,j)) > tol
                I(k,1) = j;
                u = invB * A(:,j);
                invB(k,:) = 1 / u(k,1) * invB(k,:);
                for i = 1:S(1,1)
                    if i == k
                        continue;
                    else
                        invB(i,:) = invB(i,:) - u(i,1) * invB(k,:);
                    end
                end
                iter = iter + 1;
                break;
            end
        end
    end
end
end