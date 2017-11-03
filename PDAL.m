function [xcur, ycur, epoch, fx] = PDAL(A, b, iters, acc,init)
%% PDAL from paper
% 'Malitsky, Pock - A first-order primal-dual algorithm with
% linesearch(2016)'
% implementation of algorithm 3
% solves 1/2(x,A'*A*x)-A'*b*x
% to compare, the input of corresponding CBCD functions 
% should be A'*A and A'*b
% input          A, b 
% output         x*, y(fval)
%%
% init
disp('To compare, the input of corresponding CBCD functions should be A^T*A and A^T*b');
dim  = size(A,1);
xcur = ones(dim,1)*init;
ycur = A*xcur-b;
miu  = 0.5;
tau  = 1/max(svds(A));
beta = 1;
theta= 1;
gamma= 0.5;
residual = 1;
A2 = A'* A;
b2 = A'* b;
epoch = 1;
fx=zeros(iters,1)/0;
% main iteration
while residual>acc && epoch<=iters
    % 1 compute
    x_old = xcur;
    xcur = max(min(xcur - tau*A*ycur,1),0);
    beta = beta/(1+gamma*beta*tau);
    % 2 kinesearch
    tau_old= tau;
    tau = tau*((sqrt(1+theta)-1)*rand()+1);
    % 2.a compute
    theta = tau/tau_old;
    sigma = beta*tau;
    xbar = xcur + theta*(xcur - x_old);
    y_old = ycur;
    ycur = (ycur+sigma*A*xbar-sigma*b)/(1+sigma);
    % 2.b break line if
    break_line = sqrt(beta)*tau*norm(A*(ycur-y_old))/norm(ycur-y_old);
    %fprintf('break_line = %.8f\n',break_line);
    count=1;
    while break_line>1 && count<=50
        % 2.a compute
        tau = tau*miu;
        theta = tau/tau_old;
        sigma = beta*tau;
        xbar = xcur + theta*(xcur - x_old);
        y_old = ycur;
        ycur = (ycur+sigma*A*xbar-sigma*b)/(1+sigma);
        % 2.b break line
        break_line = sqrt(beta)*tau*norm(A*(ycur-y_old))/norm(ycur-y_old);
        %fprintf('break_line = %.8f\n',break_line);
        count=count+1;
    end
    % show residual
    grad = A2*xcur;
    index_l = find(xcur<=0+2*eps);
    index_u = find(xcur>=1-2*eps);
    index = find(xcur>0+2*eps & xcur<1-2*eps);
    residual = norm([grad(index)-b2(index);min(0,grad(index_l)-b2(index_l));max(0,grad(index_u)-b2(index_u))],2);
    %yout = fval(A2,b2,xcur);
    fx(epoch)=fval(A2,b2,xcur);
    if mod(epoch, 100)==0
        fprintf('epoch;%5d, residual:%.15f, fval:%.15f\n',epoch,residual,fx(epoch));
    end
    epoch = epoch+1;
end
fx(isnan(fx))=[];
end