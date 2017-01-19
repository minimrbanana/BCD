function p=time_evaluation(d,iters,mode)
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

switch mode
    case 1
        % A tri-diagonal
        % in sparse coding
        e = ones(d,1);
        A = spdiags([-e,-e,2*e,-e,-e],[-d+1, -1, 0, 1, d-1],d,d);
    
    case 2
        % A laplacian of grid
        % matrix D
        n=sqrt(d);
        assert(n==round(n),'d is not a square number!');
        D = ones(n,n)*4;
        D(1,:)=D(1,:)-1;
        D(n,:)=D(n,:)-1;
        D(:,1)=D(:,1)-1;
        D(:,n)=D(:,n)-1;
        D = reshape(D',n*n,1);
        % matrix A
        A_below = ones(n,n);
        A_above = ones(n,n);
        A_below(:,n)=0;
        A_above(:,1)=0;
        A_below = reshape(A_below',n*n,1);
        A_above = reshape(A_above',n*n,1);
        A1 = ones(n*n,1);
        A = spdiags([-A1,-A_below,D,-A_above,-A1],[-n,-1,0,1,n],d,d);
    otherwise
        % default tri-diagonal
        e = ones(d,1);
        A = spdiags([-e,-e,2*e,-e,-e],[-d+1, -1, 0, 1, d-1],d,d);
end

% b Gaussian dist.
b = randn(d,1);

%A = sparse(A);
profile on;
% functions to evaluate
[x1sp, y1sp] = CBCD_size1_mex_sparse(A, b, d, iters);
[x2sp, y2sp] = CBCD_size2_9_mex_sparse(A, b, d, iters);
[x2sp, y2sp] = CBCD_size2_ss(A, b, d, iters);
[x3sp, y3sp] = CBCD_size3_mex_27_sparse(A, b, d, iters);
[x3sp, y3sp] = CBCD_size3_ss(A, b, d, iters);
p=profile('info');
profile off;
epochs = ones(5,1);
epochs(1) = length(y1sp);epochs(2) = length(y2sp);epochs(3) = length(y2sp);
epochs(4) = length(y3sp);epochs(5) = length(y3sp);
if size(p.FunctionTable,1)==5
    for i=1:size(p.FunctionTable,1)
        fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
        fprintf('Runtime : %.4f seconds.   ',p.FunctionTable(i,1).TotalTime);
        fprintf('T/epoch : %.4f seconds\n',p.FunctionTable(i,1).TotalTime/epochs(i));
    end
end
    
    
    
    
    