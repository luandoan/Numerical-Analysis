 function intf = simpson(f, a, b, N)

% function simpson(f, a, b, N)
% Composite Simpson's rule to approximate the integral of f from a to b.
% Uses N+1 function evaluations (N even, N>=2;  h = (b-a)/N).  

 if mod(N,2)~=0, fprintf('simpson error: N must be an even integer!\n');end

 h = (b-a)/N; intf = 0;
 intf = feval(f,a)+4*feval(f,a+h)+feval(f,b);
 for j=2:2:N-2
    intf = intf + 2*feval(f,a+j*h) + 4*feval(f,a+(j+1)*h);
 end
 intf = intf*h/3;
