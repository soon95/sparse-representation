function y = l0norm(x)
%        y = l0norm(x)
%
%  --- l0 norm ox x
y = sum(abs(x(:))~=0);
