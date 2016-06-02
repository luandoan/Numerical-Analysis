function xstar = halley(f,fprime,f2prime, x0)
 maxit = 60; 
 fx = feval(f,x0); x=x0; k=0;       % initialize
 fprintf(' %3d  %20.14f  %10.7e\n', k, x, fx);
 while (abs(fx) > 1e-15) && (k < maxit)
    x = x - (2*fx*feval(fprime,x))/(2*feval(fprime,x)*feval(fprime,x) - fx*feval(f2prime,x));
    k = k+1;  
    fx = feval(f,x);
    fprintf(' %3d  %20.14f  %10.7e\n', k, x, fx);
 end
 xstar = x;