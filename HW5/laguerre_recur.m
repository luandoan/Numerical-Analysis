function [ p2, dp2, p1 ] = laguerre_recur ( x, n, a, b, c )
  p1 = 1.0;
  dp1 = 0.0;
  p2 = x - a - 1.0;
  dp2 = 1.0;
  for i = 2 : n
        p0 = p1;    dp0 = dp1;
        p1 = p2;    dp1 = dp2;
        p2 = ( x - b(i) ) * p1 - c(i) * p0;
        dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0;
  end
end