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
diag_A1=ones(1,d-1);
A = eye(d)*2-diag(diag_A1,1)-diag(diag_A1',-1);
A(1,d)=-1;A(d,1)=-1;
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
    
    
    
    
    