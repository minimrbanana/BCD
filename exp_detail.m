function EXP = exp_detail(exp_idx)
%
% 
%

% EXP.xxx to save the parameters of the experiment

EXP.idx = exp_idx;
% default paras
EXP.isplot = 0;  % plot matrix A 
EXP.plot_convergence = 0; % plot the averaged convergence
EXP.save = 0;
EXP.max_iter = 3000;  % max number of iters
% the parameter details of each exam
rng(1);
switch exp_idx
    case 1
        % A has block size2 on the diagonal
        d = 3000;
        e0 = [1;0];
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = spdiags([-e1,-[0;e1(1:end-1)]],[-1,1],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
    case 2
        % A has block size2 on the diagonal with shift
        d = 3000;
        e0 = [0;1];
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = spdiags([-e1,-[0;e1(1:end-1)]],[-1,1],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.A(1,1)=1;EXP.A(d,d)=1;
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 3
        % A has block size2 on the diagonal
        % with noise off the diagonal
        d = 3000;
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
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 4
        % A has block size2 on the diagonal with shift
        % with noise off the diagonal
        d = 3000;
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
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 5
        % A has block size3 on the diagonal
        d = 3000;
        e0 = [1;1;0];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;0;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = spdiags([-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],[-2,-1,1,2],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 6
        % A has block size3 on the diagonal with shift
        d = 3000;
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
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;   
    case 7
        % A has block size3 on the diagonal
        % with noise off the diagonal
        d = 3000;
        e0 = [1;1;0];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;0;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e2,-e1,e1*0,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],[-2,-1,0,1,2],A_off);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 8
        % A has block size3 on the diagonal with shift
        % with noise off the diagonal
        d = 3000;
        e0 = [0;1;1];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [0;1;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e2,-e1,e1*0,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],[-2,-1,0,1,2],A_off);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 9
        % A has block size4 on the diagonal
        d = 3000;
        e0 = [1;1;1;0];
        e1 = e0(:,ones(ceil(d/4),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;1;0;0];
        e2 = e0(:,ones(ceil(d/4),1));
        e2 = reshape(e2 ,numel(e2),1);
        e0 = [1;0;0;0];
        e3 = e0(:,ones(ceil(d/4),1));
        e3 = reshape(e3 ,numel(e3),1);
        A = spdiags([-e3,-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)],-[0;0;0;e3(1:end-3)]],[-3,-2,-1,1,2,3],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 10
        % A has block size5 on the diagonal
        d = 3000;
        e0 = [0;1;1;1];
        e1 = e0(:,ones(ceil(d/4),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [0;1;1;0];
        e2 = e0(:,ones(ceil(d/4),1));
        e2 = reshape(e2 ,numel(e2),1);
        e0 = [0;1;0;0];
        e3 = e0(:,ones(ceil(d/4),1));
        e3 = reshape(e3 ,numel(e3),1);
        A = spdiags([-e3,-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)],-[0;0;0;e3(1:end-3)]],[-3,-2,-1,1,2,3],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
        EXP.save = 1;
    case 100
        % A is a tri-diagonal matrix without noise
        d = 3000; 
        e1 = ones(d,1);
        EXP.A = spdiags([-e1,-e1],[-1,1],d,d);
        diagonal = -sum(EXP.A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,EXP.A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 101
        % A is a 5-band matrix without noise
        d = 3000; 
        e1 = ones(d,1);
        EXP.A = spdiags([-e1,-e1,-e1,-e1],[-2,-1,1,2],d,d);
        diagonal = -sum(EXP.A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,EXP.A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 102
        % A is a 7-band matrix without noise
        d = 3000; 
        e1 = ones(d,1);
        EXP.A = spdiags([-e1,-e1,-e1,-e1,-e1,-e1],[-3,-2,-1,1,2,3],d,d);
        diagonal = -sum(EXP.A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,EXP.A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 200
        % A is a random symmetric matrix 
        d = 30; 
        A = sprandsym(d,3/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        EXP.A = spdiags(diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 1000;
        EXP.isplot = 0;
        EXP.plot_convergence = 0;
    case 900
        % A has only one off-entry 
        % which means after reordering, there is only one size2 block
        d = 10;
        EXP.A = speye(d);
        EXP.A(floor(d/3),floor(2*d/3))=-1;EXP.A(floor(2*d/3),floor(d/3))=-1;
        EXP.A(floor(2*d/3),floor(2*d/3))=2;EXP.A(floor(d/3),floor(d/3))=2;
        EXP.d = d;
        EXP.n_loop = 100;
        EXP.isplot = 1;
        EXP.plot_convergence = 1;
end
EXP.output_dir = ['./result/EXP_idx_' num2str(EXP.idx) '/'];
if ~exist(EXP.output_dir,'dir')
    mkdir(EXP.output_dir);
end
end