function [x,y] = RBCD_size1(A, b, d, lower, upper, max_iter)
% Random Block Coordinate Descent method to solve
% min 1/2<x,Ax>-<b,x>
% s.t. x in R^d, lower(i)<=x(i)<=upper(i)
% with block size 1
% input: A: in R^d*d
%        b: in R^d
%        d: dimension
%        lower, upper: bounds
%        max_iter: max iteration
% output: solution x

% init x and Lipschitz constant
x = lower;
y = zeros(max_iter,1);
residual = zeros(max_iter,1);
L = diag(A); % for quadratic functions the Lipschitz constant is A_ii
% pL = L;
% for i=2:d
%     pL(i) = pL(i-1)+L(i);
% end
% pL=pL/sum(L);
nL = L/sum(L);
for k=1:max_iter
    % how t choose i
    i = randsample(d,1,true,nL);
    x(i) = x(i) - (A(i,:)*x-b(i))/L(i);
    x(i) = max(lower(i),min(upper(i),x(i)));% bounds
    residual(k) = norm(A*x-b,2);
    y(k) = fval(A,b,x);
    fprintf('iter;%5d, residual:%.8f, fval:%.8f\n',k,residual(k),y(k));
end

end