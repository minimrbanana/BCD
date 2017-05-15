function EXP = exp_detail(exp_idx)
%
% 
%

% EXP.xxx to save the parameters of the experiment

EXP.idx = exp_idx;
% default paras
EXP.isplot = 0;  % plot matrix A,B & C 
EXP.plot_convergence = 0; % plot the averaged convergence
EXP.save = 0;
% max number of iters, should not set too large, 
% otherwise cannot save all the convergence matrices (dim=n_loop*mex_iter)
EXP.max_iter = 200000;  
% precision
EXP.precision = 1E-10;
% the bounds and initial state and alpha in RBCD
EXP.lower = 0;
EXP.upper = 1;
EXP.init  = 0;
EXP.alpha = 1.0;
% the parameter details of each exam
rng(1);
switch exp_idx
    case 1
        % A has block size2 on the diagonal
        d = 5000;
        e0 = [1;0];
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = spdiags([-e1,-[0;e1(1:end-1)]],[-1,1],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 2
        % A has block size2 on the diagonal with shift
        d = 5000;
        e0 = [0;1];
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = spdiags([-e1,-[0;e1(1:end-1)]],[-1,1],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.A(1,1)=1;EXP.A(d,d)=1;
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 3
        % A has block size2 on the diagonal
        % with noise off the diagonal
        d = 5000;
        e0 = [1;0];
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,e1*0,-[0;e1(1:end-1)]],[-1,0,1],A_off);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 4
        % A has block size2 on the diagonal with shift
        % with noise off the diagonal
        d = 5000;
        e0 = [0;1]; % change from [0,1] to [1,0] as shift
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,e1*0,-[0;e1(1:end-1)]],[-1,0,1],A_off);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 5
        % A has block size3 on the diagonal
        d = 5000;
        e0 = [1;1;0];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;0;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = spdiags([-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],...
            [-2,-1,1,2],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 6
        % A has block size3 on the diagonal with shift
        d = 5000;
        e0 = [0;1;1];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [0;1;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = spdiags([-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],[-2,-1,1,2],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1; 
    case 7
        % A has block size3 on the diagonal
        % with noise off the diagonal
        d = 5000;
        e0 = [1;1;0];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;0;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e2,-e1,e1*0,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],...
            [-2,-1,0,1,2],A_off);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 8
        % A has block size3 on the diagonal with shift
        % with noise off the diagonal
        d = 5000;
        e0 = [0;1;1];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [0;1;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e2,-e1,e1*0,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],...
            [-2,-1,0,1,2],A_off);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 9
        % A has block size4 on the diagonal
        d = 5000;
        e0 = [1;1;1;0];
        e1 = e0(:,ones(ceil(d/4),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;1;0;0];
        e2 = e0(:,ones(ceil(d/4),1));
        e2 = reshape(e2 ,numel(e2),1);
        e0 = [1;0;0;0];
        e3 = e0(:,ones(ceil(d/4),1));
        e3 = reshape(e3 ,numel(e3),1);
        A = spdiags([-e3,-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)],...
            -[0;0;0;e3(1:end-3)]],[-3,-2,-1,1,2,3],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 10
        % A has block size5 on the diagonal
        d = 5000;
        e0 = [0;1;1;1];
        e1 = e0(:,ones(ceil(d/4),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [0;1;1;0];
        e2 = e0(:,ones(ceil(d/4),1));
        e2 = reshape(e2 ,numel(e2),1);
        e0 = [0;1;0;0];
        e3 = e0(:,ones(ceil(d/4),1));
        e3 = reshape(e3 ,numel(e3),1);
        A = spdiags([-e3,-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)],...
            -[0;0;0;e3(1:end-3)]],[-3,-2,-1,1,2,3],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 100
        % A is a tri-diagonal matrix without noise
        d = 5000; 
        e1 = ones(d,1);
        EXP.A = spdiags([-e1,-e1],[-1,1],d,d);
        diagonal = -sum(EXP.A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,EXP.A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 101
        % A is a 5-band matrix without noise
        d = 5000; 
        e1 = ones(d,1);
        EXP.A = spdiags([-e1,-e1,-e1,-e1],[-2,-1,1,2],d,d);
        diagonal = -sum(EXP.A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,EXP.A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 102
        % A is a 7-band matrix without noise
        d = 5000; 
        e1 = ones(d,1);
        EXP.A = spdiags([-e1,-e1,-e1,-e1,-e1,-e1],[-3,-2,-1,1,2,3],d,d);
        diagonal = -sum(EXP.A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,EXP.A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 103
        % A is a tri-diagonal matrix WITH noise
        d = 5000; 
        e1 = ones(d,1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,e1*0,-e1],[-1,0,1],A_off);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 200
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,3/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 201
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,5/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 202
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,7/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 203
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,9/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 204
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,11/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;    
    case 205
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,13/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1; 
    case 206
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,15/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;  
    case 207
        % A is a random symmetric matrix 
        d = 5000; 
        A = sprandsym(d,17/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        % not clear why the precision cannot reach 1E-13
        EXP.d = d;
        EXP.n_loop = 1;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;         
    case 900
        % A has only one off-entry 
        % which means after reordering, there is only one size2 block
        % b need to be set to special value
        d = 12;
        EXP.A = speye(d);
        EXP.A(1,2)=-1;EXP.A(2,1)=-1;
        EXP.A(1,1)=2;EXP.A(2,2)=2;
        EXP.d = d;
        EXP.n_loop = 1;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
    case 901
        % A has only one off-entry and A \in 2*2
        % which means A is a size2 block
        % b need to be set to special value
        d = 4;
        EXP.A = sparse([1.1,-1,0,0;-1,1.1,0,0;0,0,1.1,-1;0,0,-1,1.1]);
        EXP.d = d;
        EXP.n_loop = 1;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;   
    case 902
        % worst case paper Ac
        d = 100;
        c = 0.8;
        A = ones(d,d)*c;
        A = A + diag(ones(d,1)*(1-c));
        EXP.A = sparse(A);
        EXP.d = d;
        EXP.n_loop = 100;
        % for the worst case the interval must be larger.
        % here we take [-10,10]
        EXP.lower = -10;
        EXP.upper = 10;
        EXP.precision = 1E-5;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
    case 903
        % worst case paper Ac with some rotation
        d = 100;
        c = 0.8;
        A = ones(d,d)*c;
        A = A + diag(ones(d,1)*(1-c));
        [U,~,~] = svd(1-2*rand(d,d));
        A = U'*A*U;
        EXP.A = sparse(A);
        EXP.d = d;
        EXP.n_loop = 100;
        % for the worst case the interval must be larger.
        % here we take [-10,10]
        EXP.lower = -10;
        EXP.upper = 10;
        EXP.precision = 1E-5;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
        EXP.save = 1;
end
EXP.output_dir = ['./result/EXP_idx_' num2str(EXP.idx) '/'];
if ~exist(EXP.output_dir,'dir')
    mkdir(EXP.output_dir);
end
end