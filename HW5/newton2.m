 function xstar = newton2(f,fprime,f2prime, x0)

% function xstar = newton(f,fprime,x0)
% Compute a root of the function f using Newton's method
% f:      a function name
% fprime: a derivative function name
% x0:     the starting guess
% Example: newton('sin','cos',3), or newton('my_f','my_fprime',1)

 maxit = 60; 
 fx = feval(f,x0); x=x0; k=0;       % initialize
 fprintf(' %3d  %20.14f  %10.7e\n', k, x, fx);
 while (abs(fx) > 1e-15) & (k < maxit)
    a = feval(f2prime,x);
    b = 2*feval(f2prime,x)*x - 2*feval(fprime,x);
    c = 2*fx - 2*feval(fprime,x)*x;
    delta = b.^2 - 4*a*c;
    x1 = (b-sqrt(delta))/a;
    x2 = (b+sqrt(delta))/a;
    
    x = x2; % - fx/feval(fprime,x);     % Newton's method
    k = k+1;  
    fx = feval(f,x);
    fprintf(' %3d  %20.14f  %10.7e\n', k, x, fx);
 end
 xstar = x;
