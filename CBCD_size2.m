function [x,y] = CBCD_size2(A, b, d, lower, upper, max_iter)
% Cyclic Block Coordinate Descent method to solve
% min 1/2<x,Ax>-<b,x>
% s.t. x in R^d, lower(i)<=x(i)<=upper(i)
% with block size 2
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
residual = zeros(max_iter,1);
index_0 = find(x==0);
index_1 = find(x==1);
res_vec = -A*x+b;
res_vec(index_0) = min(0,A(index_0,:)*x+b(index_0));
res_vec(index_1) = max(0,A(index_1,:)*x+b(index_1));
residual(1) = norm(res_vec,2);
fprintf('epoch;    0, residual:%.15f, fval:%.8f\n',residual(1),y(1));
epoch = 1;
while residual(epoch)~=0 && epoch<=max_iter
    for i=1:floor(d/2)
        % formulate size 2 problem
        % index of loop
        index = i*2-1;
        % A2 is the index and index+1 th row and column of A
        A2 = [A(index,index),A(index,index+1);A(index+1,index),A(index+1,index+1)];
        x2 = x(index:index+1,1);
        % b2 is b, minus inner product of A(index(or +1),:),x and
        % A(:,index(or +1)),x
        % !!!without term in A2!!!
        b2 = [b(index)-(A(index,:)*x/2 - A2(1,:)*x2)-(A(:,index)'*x/2 - A2(:,1)'*x2);b(index+1)-(A(index+1,:)*x/2 - A2(2,:)*x2)-(A(:,index+1)'*x/2 - A2(:,2)'*x2)];
        % 9 cases to choose
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

