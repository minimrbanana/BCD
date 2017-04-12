function [x,residual] = RBCD_size1(A, b, d, lower, upper, max_iter)
% Random Block Coordinate Descent method to solve
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
residual = ones(max_iter,1);
grad = A*x;
index_l = find(x<=lower+2*eps);
index_u = find(x>=upper-2*eps);
index = find(x>lower+2*eps & x<upper-2*eps);
residual(1) = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));max(0,grad(index_u)-b(index_u))],2);
fprintf('epoch;    0, residual:%.15f, fval:%.8f\n',residual(1),y(1));
epoch = 1;
nL = L/sum(L);
% set random seed
stream = RandStream.getGlobalStream;
reset(stream);
while residual(epoch)>1E-13 && epoch<=max_iter
    for k=1:d
        % how t choose i
        i = randsample(d,1,true,nL);
        %fprintf('i;%5d\n',i);
        x(i) = x(i) - (A(i,:)*x-b(i))/L(i);
        x(i) = max(lower(i),min(upper(i),x(i)));% bounds
    end
    % compute the real gradient after each epoch
    grad = A*x;
    % opt condition, 0 in sub gradient
    index_l = find(x<=lower+2*eps);
    index_u = find(x>=upper-2*eps);
    index = find(x>lower+2*eps & x<upper-2*eps);
    residual(epoch+1) = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));max(0,grad(index_u)-b(index_u))],2);
    if(rem(epoch,4)==0)
        fprintf('epoch;%5d, residual:%.15f\n',epoch,residual(epoch+1));
    end
    epoch = epoch+1;
end
% show residual of last epoch
fprintf('epoch;%5d, residual:%.15f\n',epoch-1,residual(epoch));
% output, cut the unvalued residual
residual(epoch+1:end)=[];

end