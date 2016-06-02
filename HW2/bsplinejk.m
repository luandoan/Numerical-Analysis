 function Bjkx = bsplinejk(x, j, k, knots, x0indx)

% function Bjkx = bsplinejk(x, j, k, knots, x0indx);
%
% B-Spline evaluation routine:  computes B_j^k(x).
%
% x: value at which to evaluate B_j^k; x can be a vector.
% j: the spline is supported on (x_j,x_{j+k+1}) for k>0, and on [x_j,x_{j+1}) for k=0.
% k: degree of polynomials that make up the spline.
% knots:  set of knot values.  
%         Must be large enough to contain the support of B_j^k.
% x0indx: array index indicating where x0 is located in knots:  knots(x0indx) = x0.

 if k==0                                                                                                        
   Bjkx = (x>=knots(x0indx+j))&(x<knots(x0indx+j+1));
 else
   Bjkx = (x-knots(x0indx+j))/(knots(x0indx+j+k)-knots(x0indx+j)).*Bx(x,j,k-1,knots,x0indx) + ...
          (knots(x0indx+j+k+1)-x)/(knots(x0indx+j+k+1)-knots(x0indx+j+1)).*Bx(x,j+1,k-1,knots,x0indx);
 end 
