function G = gauss_laguerre(n,z)
a = z-1;
for i =1:n
    b(i) = (a+2*i-1);
    c(i) = (i-1)*(a+i-1);
end
cc = gamma(a+1)*prod(c(2:n));
for i = 1:n
    [x,dp2,p1 ] = laguerre_root(x,n,a,b,c);
    xt(i) = x;
    w(i) = (cc/dp2)/p1;
    if (i==1) 
        x = ((1+a)*(3+0.92*a))/(1+2.4*n+1.8*a);
    elseif (i==2)
        x = x + (15+6.25*a)/(1+0.9*a+2.5*n);
    else
        r1 = (1 + 2.55*(i-2))/(1.9*(i-2));
        r2 = 1.26 *(i-2)*a /(1+3.5*(i-2));
        ratio = (r1+r2)/(1+0.3*a);
        x = x + ratio*(x-xt(i-2) );
    end
    
%     [x,dp2,p1 ] = laguerre_root(x,n,a,b,c);
%     xt(i) = x;
%     w(i) = (cc/dp2)/p1;
end
  G = 0;
  for i = 1:n
      G = G + w(i)*xt(i).^a;
  end
end