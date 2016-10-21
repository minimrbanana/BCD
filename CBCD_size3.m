function [x,y] = CBCD_size3(A, b, dim, lower, upper, max_iter)
% Cyclic Block Coordinate Descent method to solve
% min 1/2<x,Ax>-<b,x>
% s.t. x in R^d, lower(i)<=x(i)<=upper(i)
% with block size 3
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
% for computing residual, based on the normal cone
residual = ones(max_iter+1,1);
fprintf('epoch;    0, residual:%.15f, fval:%.15f\n',residual(1),y(1));
epoch = 1;
grad = A*x;
while residual(epoch)>1E-13 && epoch<=max_iter
    for ii=1:floor(dim/3)
        i = ii*3-2;
        %grad = grad-A(:,[i,i+1])*x([i,i+1]);
        grad = grad-A(:,i)*x(i)-A(:,i+1)*x(i+1)-A(:,i+2)*x(i+2);
        % update x(i)
        % define size 2 block
        %A2 = [A(i,i),A(i,i+1),A(i,i+2);A(i+1,i),A(i+1,i+1),A(i+1,i+2);A(i+2,i),A(i+2,i+1),A(i+2,i+2)];
        %b2 = [b(i);b(i+1);b(i+2)]-[grad(i);grad(i+1);grad(i+2)];
        a11 = A(i,i);   a12 = A(i,i+1);   a13 = A(i,i+2);
        a21 = A(i+1,i); a22 = A(i+1,i+1); a23 = A(i+1,i+2);
        a31 = A(i+2,i); a32 = A(i+2,i+1); a33 = A(i+2,i+2);
        b1 = b(i)-grad(i); b2 = b(i+1)-grad(i+1); b3 = b(i+2)-grad(i+2);
        % decission tree begin
        flag = 0;
        % first discuss three 0 or 1, 8 cases
        if b1<=0 && b2<=0 && b3<=0 % case1
            x([i,i+1,i+2])=[0;0;0];flag = 1;
        else if a13>=b1 && a23>=b2 && a33<=b3 % case3
                x([i,i+1,i+2])=[0;0;1];flag = 1;
            else if a12>=b1 && a22<=b2 && a32>=b3 % case7
                    x([i,i+1,i+2])=[0;1;0];flag = 1;
                else if a12+a13>=b1 && a22+a23<=b2 && a32+a33<=b3 % case9
                        x([i,i+1,i+2])=[0;1;1];flag = 1;
                    else if a11<=b1 && a12>=b2 && a13>=b3 % case19
                            x([i,i+1,i+2])=[1;0;0];flag = 1;
                        else if a11+a13<=b1 && a21+a23>=b2 && a31+a33<=b3 % case21
                                x([i,i+1,i+2])=[1;0;1];flag = 1;
                            else if a11+a12<=b1 && a21+a22<=b2 && a31+a32>=b3 % case25
                                    x([i,i+1,i+2])=[1;1;0];flag = 1;
                                else if a11+a12+a13<=b1 && a21+a22+a23<=b2 && a31+a32+a33<=b3 % case27
                                        x([i,i+1,i+2])=[1;1;1];flag = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        % second discuss two 0s and two 1s, 6 cases
        if flag==0
            x3 = b3/a33;
            if x3>=0 && x3<=1 && a13*x3>=b1 && a23*x3>=b2 % case2
                x([i,i+1,i+2])=[0;0;x3];flag = 1;
            end
            if flag==0
                x2 = b2/a22;
                if x2>=0 && x2<=1 && a12*x2>=b1 && a32*x2>=b3 % case4
                    x([i,i+1,i+2])=[0;x2;0];flag = 1;
                end
                if flag==0
                    x1 = b1/a11;
                    if x1>=0 && x1<=1 && a21*x1>=b2 && a31*x1>=b3 % case10
                        x([i,i+1,i+2])=[x1;0;0];flag = 1;
                    end
                    if flag==0
                        x1=(b1-a12-a13)/a11;
                        if x1>=0 && x1<=1 && a21*x1+a22+a23<=b2 && a31*x1+a32+a33<=b3 % case18
                            x([i,i+1,i+2])=[x1;1;1];flag = 1;
                        end
                        if flag==0
                            x2=(b2-a21-a23)/a22;
                            if x2>=0 && x2<=1 && a11+a12*x2+a13<=b1 && a31+a32*x2+a33<=b3 % case24
                                x([i,i+1,i+2])=[1;x2;1];flag = 1;
                            end
                            if flag==0
                                x3=(b3-a31-a32)/a33;
                                if x3>=0 && x3<=1 && a11+a12+a13*x3<=b1 && a21+a22+a23*x3<=b2 % case26
                                   x([i,i+1,i+2])=[1;1;x3];flag = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        % third discuss one 0 and one 1
        if flag==0
            x1=(b1-a13)/a11;
            if x1>=0 && x1<=1 && a21*x1+a23>=b2 && a31*x1+a33<=b3 % case12
                x([i,i+1,i+2])=[x1;0;1];flag = 1;
            end
            if flag==0
                x1 = (b1-a12)/a11;
                if x1>=0 && x1<=1 && a21*x1+a22<=b2 && a31*x1+a32>=b3 % case16
                    x([i,i+1,i+2])=[x1;1;0];flag = 1;
                end
                if flag==0
                    x2 = (b2-a23)/a22;
                    if x2>=0 && x2<=1 && a12*x2+a13>=b1 && a32*x2+a33<=b3 % case6
                        x([i,i+1,i+2])=[0;x2;1];flag = 1;
                    end
                    if flag==0
                        x2=(b2-a21)/a22;
                        if x2>=0 && x2<=1 && a11+a12*x2<=b1 && a31+a32*x2>=b3 % case22
                            x([i,i+1,i+2])=[1;x2;0];flag = 1;
                        end
                        if flag==0
                            x3=(b3-a32)/a33;
                            if x3>=0 && x3<=1 && a12+a13*x3>=b1 && a22+a23*x3<=b2 % case8
                                x([i,i+1,i+2])=[0;1;x3];flag = 1;
                            end
                            if flag==0
                                x3=(b3-a31)/a33;
                                if x3>=0 && x3<=1 && a11+a13*x3<=b1 && a21+a23*x3>=b2 % case26
                                   x([i,i+1,i+2])=[1;0;x3];flag = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        % fourth discuss one 1 or one 0
        if flag==0
            x23=[a22,a23;a32,a33]\[b2;b3];
            if x23(1)>=0 && x23(1)<=1 && x23(2)>=0 && x23(2)<=1 && a12*x23(1)+a13*x23(2)>=b1 % case5
                x([i,i+1,i+2])=[0;x23];flag = 1;
            end
            if flag==0
                x23=[a22,a23;a32,a33]\[b2-a21;b3-a31];
                if x23(1)>=0 && x23(1)<=1 && x23(2)>=0 && x23(2)<=1 && a11+a12*x23(1)+a13*x23(2)<=b1 % case23
                    x([i,i+1,i+2])=[1;x23];flag = 1;
                end
                if flag==0
                    x13=[a11,a13;a31,a33]\[b1;b3];
                    if x13(1)>=0 && x13(1)<=1 && x13(2)>=0 && x13(2)<=1 && a21*x13(1)+a23*x13(2)>=b2 % case11
                        x([i,i+1,i+2])=[x13(1);0;x13(2)];flag = 1;
                    end
                    if flag==0
                        x13=[a11,a13;a31,a33]\[b1-a12;b3-a32];
                        if x13(1)>=0 && x13(1)<=1 && x13(2)>=0 && x13(2)<=1 && a21*x13(1)+a22+a23*x13(2)<=b2 % case17
                            x([i,i+1,i+2])=[x13(1);1;x13(2)];flag = 1;
                        end
                        if flag==0
                            x12=[a11,a12;a21,a22]\[b1;b2];
                            if x12(1)>=0 && x12(1)<=1 && x12(2)>=0 && x12(2)<=1 && a31*x12(1)+a32*x12(2)>=b3 % case13
                                x([i,i+1,i+2])=[x12;0];flag = 1;
                            end
                            if flag==0
                                x12=[a11,a12;a21,a22]\[b1-a13;b2-a23];
                                if x12(1)>=0 && x12(1)<=1 && x12(2)>=0 && x12(2)<=1 && a31*x12(1)+a32*x12(2)+a33<=b3 % case15
                                    x([i,i+1,i+2])=[x12;1];flag = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        % last case, no 0 or 1
        if flag==0
            x123=[a11,a12,a13;a21,a22,a23;a31,a32,a33]\[b1;b2;b3];
            if x123(1)>=0 && x123(1)<=1 && x123(2)>=0 && x123(2)<=1 && x123(3)>=0 && x123(3)<=1 % case14
                x([i,i+1,i+2])=x123;flag = 1;
            else 
                fprintf('no update, check code\n');
            end
        end
        
        %fprintf('x_%d = %.15f, x_%d = %.15f\n',i,x(i),i+1,x(i+1));
        % decission tree end
        %pause;
        grad = grad+A(:,i)*x(i)+A(:,i+1)*x(i+1)+A(:,i+2)*x(i+2);
    end
    % if mod(dim,2)==1
    if mod(dim,3)==1
        i=dim;
        grad = grad-A(:,i)*x(i);
        x(i) = max(lower(i),min(upper(i),(b(i)-grad(i))/A(i,i)));
        grad = grad+A(:,i)*x(i);
        
    else if mod(dim,3)==2
        i=dim-1;
        grad = grad-A(:,i)*x(i);
        x(i) = max(lower(i),min(upper(i),(b(i)-grad(i))/A(i,i)));
        grad = grad+A(:,i)*x(i);
        
        i=dim;
        grad = grad-A(:,i)*x(i);
        x(i) = max(lower(i),min(upper(i),(b(i)-grad(i))/A(i,i)));
        grad = grad+A(:,i)*x(i);
        
        end
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