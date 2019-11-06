function z = hard(x,lambda)
%        z = hard(x,lambda)
%
%  --- Hard threshold ----
%
%  z is the solution of 
%
%      z= arg min_x = 0.5*|| y - x ||_2^2 + lambda || x ||_0, 
%

z = x.*(abs(x) >= sqrt(2*lambda));
