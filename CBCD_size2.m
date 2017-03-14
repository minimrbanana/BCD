function [x,residual] = CBCD_size2(A, b, dim, lower, upper, max_iter)
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
%         residual r in each epoch
fprintf('CBCD size 2.m\n');
% init x and Lipschitz constant
x = lower;

% for computing residual, based on the normal cone
residual = ones(max_iter+1,1);
residual(1) = norm(min(0,-b),2);
fprintf('epoch;    0, residual:%.15f\n',residual(1));
epoch = 1;
grad = A*x;
while residual(epoch)>1E-13 && epoch<max_iter
    for ii=1:floor(dim/2)
        i = ii*2-1;
        %grad = grad-A(:,[i,i+1])*x([i,i+1]);
        grad = grad-A(:,i)*x(i)-A(:,i+1)*x(i+1);
        % update x(i)
        % define size 2 block
        A2 = [A(i,i),A(i,i+1);A(i+1,i),A(i+1,i+1)];
        b2 = [b(i);b(i+1)]-[grad(i);grad(i+1)];
        % decission tree begin
        flag = 0;
        a21b1_a11 = A2(2,1)*b2(1)/A2(1,1);% temp var1
        a12b2_a22 = A2(1,2)*b2(2)/A2(2,2);% temp var2
        if b2(1)<=0 % case b2(i)<=0
            if b2(2)<=0
                x([i,i+1])=[0;0];flag=1;%disp('case1'); % case 1
            end
        else % case b2(i)>0
            if b2(1)>=A2(1,1)
                if b2(2)<=A2(2,1)
                    x([i,i+1])=[1;0];flag=1; %disp('case3');% case 3
                end
            else
                if b2(2)<=a21b1_a11 %use temp var1
                    x([i,i+1])=[b2(1)/A2(1,1);0];flag=1; %disp('case2');% case 2
                end
            end
        end
        if flag==0 % x2~=0, assume x2=1
            if b2(1)<=A2(1,2)
                if b2(2)>=A2(2,2)
                    x([i,i+1])=[0;1];flag=1;%disp('case7');% case 7
                end
            else
                if b2(1)>=A2(1,1)+A2(1,2)
                    if b2(2)>=A2(1,2)+A2(2,2)
                        x([i,i+1])=[1;1];flag=1;%disp('case9');% case 9
                    end
                else
                    if b2(2)>=a21b1_a11+det(A2)/A2(1,1)%use temp var1
                         x([i,i+1])=[(b2(1)-A2(1,2))/A2(1,1);1];flag=1;%disp('case8'); % case 8
                    end
                end
            end
        end
        if flag==0 % x2~=0 & x2~=1, x2 in (0,1)
            if b2(1)<=a12b2_a22 %use temp var2
                x([i,i+1])=[0;b2(2)/A2(2,2)];flag=1; %disp('case4');% case 4
            else
                if b2(1)>=a12b2_a22+det(A2)/A2(2,2)%use temp var2
                    x([i,i+1])=[1;(b2(2)-A2(2,1))/A2(2,2)];flag=1; %disp('case6');% case 6
                else
                    x([i,i+1])=min(max(A2\b2,0),1);flag=1;%disp('case5'); % case 5
                end
            end
        end
        % decission tree end
        grad = grad+A(:,i)*x(i)+A(:,i+1)*x(i+1);
    end
    % if mod(dim,2)==1
    if mod(dim,2)==1
        i=dim;
        grad = grad-A(:,i)*x(i);
        x(i) = max(lower(i),min(upper(i),(b(i)-grad(i))/A(i,i)));
        % grad = grad+A(:,i)*x(i);
    end
    %compute the real gradient after each epoch
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