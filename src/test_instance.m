function [A,b,c,I] = test_instance(m,n,d)
A = zeros(m,n);
for i = 1:m
    for j = 1:n
        if rand() < d
            A(i,j) = 2 * rand() - 1;
        end
    end
end
A(:,n - m + 1:end) = eye(m);
x = rand(n,1);
b = A * x;
c = rand(n,1);
c(n - m + 1:end) = 0;
I = (n - m + 1:n)';
end