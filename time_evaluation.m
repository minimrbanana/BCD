function p=time_evaluation(d,iters)
% time test
% input : d   dimension of the QP problem
% outout: run time of functions

% set the default dim and max iterations
if nargin == 0
    d=1000; iters=200;
else if nargin == 1
        iters=200;
    end
end

% A tri-diagonal
% in sparse coding
e = ones(d,1);
A = spdiags([-e,-e,2*e,-e,-e],[-d+1 -1 0 1 d-1],d,d);


% A laplacian of grid
% matrix D
% n=sqrt(d);
% D = ones(n,n)*4;
% D(1,:)=D(1,:)-1;
% D(n,:)=D(n,:)-1;
% D(:,1)=D(:,1)-1;
% D(:,n)=D(:,n)-1;
% D = reshape(D',n*n,1);
% D = diag(D);
% % matrix A
% A = ones(n,n);
% A(:,n)=0;
% A = reshape(A',n*n,1);
% A = diag(A(1:end-1),1)+diag(A(1:end-1),-1)+ ...
%     diag(ones(n*(n-1),1),-n)+diag(ones(n*(n-1),1),n);
% 
% A=D-A;
% b Gaussian dist.
b = randn(d,1);

A = sparse(A);
profile on;
% functions to evaluate
[x1sp, y1sp] = CBCD_size1_mex_sparse(A, b, d, iters);
[x2sp, y2sp] = CBCD_size2_9_mex_sparse(A, b, d, iters);
[x3sp, y3sp] = CBCD_size3_mex_27_sparse(A, b, d, iters);
p=profile('info');
profile off;
for i=1:size(p.FunctionTable,1)
    fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
    fprintf('Runtime : %.4f seconds\n',p.FunctionTable(i,1).TotalTime);
end
    
    
    
    
    