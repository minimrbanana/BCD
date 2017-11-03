function thesis_plot(idx)
%% thesis plot
close all;
lambda=1E-5;
d = 30;
Fsize = 23;
saveDir = ['./thesis_plot/EXP_idx_' num2str(idx) '/'];
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

switch idx
    case 1
        e0 = [1;1;0];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;0;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = spdiags([-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],...
            [-2,-1,1,2],d,d);
        diagonal = -sum(A);% in order to make matrix A positive
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal'+lambda,0,A);
    case 2
        e1 = ones(d,1);
        A = spdiags([-e1,-e1],[-1,1],d,d);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal'+lambda,0,A);
    case 3
        A = sprandsym(d,3/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal'+lambda,0,A);
    case 4
        A = sprandsym(d,9/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A);
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal'+lambda,0,A);
    case 5
        c = 0.8;
        A = ones(d,d)*c;
        A = A + diag(ones(d,1)*(1-c));
        A = sparse(A);
    case 6
        c = 0.8;
        A = ones(d,d)*c;
        A = A + diag(ones(d,1)*(1-c));
        [Q,~,~] = svd(rand(d,d));
        A = Q'*A*Q;
        A = sparse(A);
        %disp('not complete');
    otherwise
        error('index not defined');
        
end

perm = randperm(d);
B = A(perm,perm);
reod = symrcm(B);
C = B(reod,reod);

figure(1),clf;
h1 = imagesc(A);
colorbar; hold on;
title('Matrix A');
set(gca,'fontsize',Fsize,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];  
saveas(gcf,[saveDir 'figure1.pdf']);

figure(2),clf;
h2 = imagesc(B);
colorbar; hold on;
title('Matrix B');
set(gca,'fontsize',Fsize,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];  
saveas(gcf,[saveDir 'figure2.pdf']);

figure(3),clf;
h3 = imagesc(C);
colorbar; hold on;
title('Matrix C');
set(gca,'fontsize',Fsize,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];  
saveas(gcf,[saveDir 'figure3.pdf']);
