function [x,V,w] = gauss_lag(n)
J = zeros(n+1);
x = zeros(1,n+1);
for i = 1:n+1
    for j = 1:n+1
        if i-j==0
            J(i,j) = 2*i-1;
        elseif j-i==1
            J(i,j) = i;
        elseif i-j==1
            J(i,j) = j;
        end
    end
end
[V,D] = eig(J);
for i = 1:size(D,1)
    for j = 1:size(D,2)
        if i-j==0
            x(1,i) = D(i,j);
        end
    end
end

w = V(:,1).^2;
