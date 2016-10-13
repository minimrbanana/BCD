function [x,y] = RBCD_size1(A, b, d, lower, upper, max_iter)
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
residual = zeros(max_iter,1);
index_0 = find(x==0);
index_1 = find(x==1);
res_vec = -A*x+b;
res_vec(index_0) = min(0,A(index_0,:)*x+b(index_0));
res_vec(index_1) = max(0,A(index_1,:)*x+b(index_1));
residual(1) = norm(res_vec,2);
fprintf('epoch;    0, residual:%.15f, fval:%.8f\n',residual(1),y(1));
epoch = 1;
nL = L/sum(L);
% set random seed
stream = RandStream.getGlobalStream;
reset(stream);
while residual(epoch)~=0 && epoch<=max_iter
    for k=1:d
        % how t choose i
        i = randsample(d,1,true,nL);
        %fprintf('i;%5d\n',i);
        x(i) = x(i) - (A(i,:)*x-b(i))/L(i);
        x(i) = max(lower(i),min(upper(i),x(i)));% bounds
    end
    % opt condition, 0 in sub gradient
    index_0 = find(x==0);
    index_1 = find(x==1);
    res_vec = -A*x+b;
    res_vec(index_0) = min(0,A(index_0,:)*x-b(index_0));
    res_vec(index_1) = max(0,A(index_1,:)*x-b(index_1));
    residual(epoch+1) = norm(res_vec,2);
    y(epoch+1) = fval(A,b,x);
    fprintf('epoch;%5d, residual:%.15f, fval:%.8f\n',epoch,residual(epoch+1),y(epoch+1));
    epoch = epoch+1;
end

end