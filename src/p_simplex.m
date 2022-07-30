function [X,optvalue,tpresolve,tphaseI,tphaseII,iter1,iter2] = p_simplex(A,C,b)
A = full(A);
S = size(A);
X = zeros(S(1,2),1);
fprintf("Presolving...\n");
tic;
[A1,b1,Ix,valx,flag] = presolve(A,b);
tpresolve = toc;
clearvars A;
if flag == 1
    return;
end
X(Ix,1) = valx;
fprintf("Finish presolving!\n");
fprintf("Starting to find the initial BFS...\n");
tic;
[I,iter1] = phaseI(A1,b1);
tphaseI = toc;
fprintf("Finish finding the initial BFS, starting to find the optimal solution...\n");
flagc = 1;
for i = 1:length(C)
    if ismember(i,Ix)
        continue;
    else
        C1(flagc,1) = C(i,1);
        flagc = flagc + 1;
    end
end
tic;
[X1,optvalue,iter2] = revised_simplex(A1,b1,C1,I);
tphaseII = toc;
flagx = 1;
for i = 1:length(X)
    if ismember(i,Ix)
        optvalue = optvalue + C(i,1) * X(i,1);
        continue;
    else
        X(i,1) = X1(flagx,1);
        flagx = flagx + 1;
    end
end
X(X < 0) = 0;
fprintf("Optimal solution found!\n");
end