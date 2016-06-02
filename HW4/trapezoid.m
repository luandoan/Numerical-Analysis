 function intf = trapezoid(f, a, b, N)

% function trapezoid(f, a, b, N)
% Composite trapezoid rule to approximate the integral of f from a to b.
% Uses N+1 function evaluations (N>1; h = (b-a)/N).

 h = (b-a)/N; intf = 0;
 intf = (feval(f,a)+feval(f,b))/2;
 for j=1:N-1
    intf = intf + feval(f,a+j*h);
 end
 intf = intf*h;
