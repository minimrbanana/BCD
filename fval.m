function y = fval(A,b,x)
% calc the function value of 
% y = 1/2<x,Ax>-<b,x>
y = x'*A*x/2-b'*x;