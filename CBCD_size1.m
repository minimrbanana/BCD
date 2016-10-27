function [x,y] = CBCD_size1(A, b, dim, lower, upper, max_iter)
% Cyclic Block Coordinate Descent method to solve
% min 1/2<x,Ax>-<b,x>
% s.t. x in R^d, lower(i)<=x(i)<=upper(i)
% with block size 1
% input:  A: in R^d*d
%         b: in R^d
%         d: dimension
%         lower, upper: bounds
%         max_iter: max iteration
% output: solution x
%         function value y in each epoch

% init x and Lipschitz constant
x = lower;
y = zeros(max_iter+1,1);
y(1) = fval(A,b,x);
L = diag(A); % for quadratic functions the Lipschitz constant is A_ii
% for computing residual, based on the normal cone
residual = ones(max_iter+1,1);
fprintf('epoch;    0, residual:%.15f, fval:%.15f\n',residual(1),y(1));
epoch = 1;
grad = A*x;
while residual(epoch)>1E-13 && epoch<=max_iter
    for i=1:dim
        grad = grad-A(:,i)*x(i);
        x(i) = max(lower(i),min(upper(i),(b(i)-grad(i))/L(i)));
        grad = grad+A(:,i)*x(i);
    end
    %compute the real gradient after each epoch
    grad = A*x;
    % opt condition, 0 in sub gradient
    index_l = find(x<=lower+2*eps);
    index_u = find(x>=upper-2*eps);
    index = find(x>lower+2*eps & x<upper-2*eps);
    residual(epoch+1) = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));max(0,grad(index_u)-b(index_u))],2);
    y(epoch+1) = fval(A,b,x);
    if(rem(epoch,1)==0)
        fprintf('epoch;%5d, residual:%.15f, fval:%.15f\n',epoch,residual(epoch+1),y(epoch+1));
    end
    epoch = epoch+1;
end
y(epoch+1:end)=[];
end

