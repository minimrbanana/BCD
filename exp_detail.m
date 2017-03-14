function EXP = exp_detail(exp_idx)
%
% 
%

% EXP.xxx to save the parameters of the experiment

EXP.idx = exp_idx;
% default paras
EXP.isplot = 0;  % plot matrix A 
EXP.plot_convergence = 0; % plot the averaged convergence
EXP.max_iter = 3000;  % max number of iters
% the parameter details of each exam
rng(1);
switch exp_idx
    case 1
        % A has block size2 on the diagonal
        d = 1000;
        e0= [1;0];
        e = e0(:,ones(floor(d/2),1));
        e = reshape(e ,numel(e),1);
        A = spdiags([-e,-[0;e(1:end-1)]],[-1,1],d,d);
        diagonal = sum(A);
        EXP.A = spdiags(-diagonal',0,A);
        EXP.d = d;
        EXP.n_loop = 40;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
    case 100
        % A is a tri-diagonal matrix without noise
        d = 5000; 
        e = ones(d,1);
        EXP.A = spdiags([-e,2*e,-e],[-1, 0, 1],d,d);
        EXP.d = d;
        EXP.n_loop = 2;
        EXP.isplot = 0;
        EXP.plot_convergence = 1;
    case 101
        d = 2000;
        EXP.A = sprandsym(d,0.09,0.5,1);
        EXP.d = d;
        EXP.n_loop = 1;
        EXP.isplot = 1;
end
EXP.output_dir = ['./result/EXP_idx_' num2str(EXP.idx) '/'];
if ~exist(EXP.output_dir,'dir')
    mkdir(EXP.output_dir);
end
end