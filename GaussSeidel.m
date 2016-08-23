function x = GaussSeidel(A, b, d, lower, upper, max_iter)
% Gauss Seidel method to solve
% min 1/2<x,Ax>-<b,x>
% s.t. x in R^d, lower(i)<=x(i)<=upper(i)
% input: A: in R^d*d
%        b: in R^d
%        d: dimension
%        lower, upper: bounds
%        max_iter: max iteration
% output: solution x
x = lower;
y = zeros(max_iter,1);
residual = zeros(max_iter,1);
for k=1:max_iter
    i = mod(k-1,d)+1;
    x(i)=max(lower(i),min(upper(i),(b(i)-A(i,:)*x+A(i,i)*x(i))/A(i,i)));
    residual(k) = norm(A*x-b,2);
    y(k) = fval(A,b,x);
    fprintf('iter;%5d, residual:%.8f, fval:%.8f\n',k,residual(k),y(k));
end
figure(1),
plot(1:max_iter,y);
hold on;
end