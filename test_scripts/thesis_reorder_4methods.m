%% reorder and BFGS, PDAL, AGD
% compare the convergence for other methods before and after reordering
addpath ../

%% type II case
% A & b
rng(1);

load /home/yu/datasets/SNAP/Social_networks/mat/facebook_combined.mat; iters=42000;
%load /home/yu/datasets/SNAP/Collaboration_networks/mat/ca-CondMat.mat; iters=400;
%load /home/yu/datasets/SNAP/Social_networks/mat/gplus_combined.mat; iters=2000;
d=size(A,1);
% d=4000;
% e1=ones(d,1);
% A=spdiags([-e1,2*e1,-e1],[-1,0,1],d,d);
% p=randperm(d,1);A=A(p,p);
A=A+speye(d)*1E-4;
xmin=rand(d,1)*2-0.5;
%xmin=rand(d,1)*0.8+0.1;
b=A*xmin;
init=0.5;
acc=1E-10;
% reorder
p1=symamd(A);
A1=A(p1,p1);b1=b(p1);
p2=symrcm(A);
A2=A(p2,p2);b2=b(p2);
% solve 


f_pdal=0;
t6=0;

%% CBCD
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,acc,0,1,init);
t1=toc(t0);

%% LBFGS-B
addpath /home/yu/LBSGF
addpath /home/yu/LBSGF/L-BFGS-B-C/
addpath /home/yu/LBSGF/L-BFGS-B-C/Matlab/
addpath /home/yu/LBSGF/L-BFGS-B-C/src/
disp('=== My test function, 2D === ');
n = d;
f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;
% constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;
% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
%trueSoln = xmin;
trueSoln = cx1;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = acc;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 20000;
% The {f,g} is another way to call it
t0=tic;
[x,f,info] = lbfgsb( {f,g} , l, u, opts );
t4=toc(t0);

if abs(f) < 1e-8
    disp('Success!');
% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
else
    disp('Something didn''t work right :-(  ');
end
fval=info.err(:,1);
f1=fval;

%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,acc,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
% modify A and b
[V,D]=eig(full(A));
A_pdal=V*spdiags(sqrt(spdiags(D)),0,d,d)*V';
b_pdal=A_pdal\b;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


fmin=0.5*trueSoln'*A*trueSoln-b'*trueSoln;
%fmin=min([cy1;fval;fx;f_pdal]);
figure(1),clf,
% semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
% semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
% semilogy(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
% semilogy(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
loglog(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
loglog(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
loglog(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
loglog(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
l1=sprintf('CBCD1,%.4f s, #%d',t1,size(cy1,1)-1);
l4=sprintf('L-BFGS-B,%.4f s, #%d',t4,size(fval,1)-1);
l5=sprintf('APGD,%.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL,%.4f s, #%d',t6,size(f_pdal,1)-1);

legend(l6,l1,l4,l5);
title('Convergence of CBCD1, L-BFGS-B, APGD and PDAL');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
axis([0 9000 1E-10 3E5]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'fontsize',16,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_reorder0.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_reorder0.mat',...
    'f_pdal','cy1','fval','fx','fmin','t1','t4','t5','t6');
%% after rcm
A=A1;b=b1;
%% CBCD
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,acc,0,1,init);
t1=toc(t0);

%% LBFGS-B
addpath /home/yu/LBSGF
addpath /home/yu/LBSGF/L-BFGS-B-C/
addpath /home/yu/LBSGF/L-BFGS-B-C/Matlab/
addpath /home/yu/LBSGF/L-BFGS-B-C/src/
disp('=== My test function, 2D === ');
n = d;
f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;
% constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;
% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
%trueSoln = xmin(p1);
trueSoln = cx1;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = acc;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 20000;
% The {f,g} is another way to call it
t0=tic;
[x,f,info] = lbfgsb( {f,g} , l, u, opts );
t4=toc(t0);

if abs(f) < 1e-8
    disp('Success!');
% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
else
    disp('Something didn''t work right :-(  ');
end
fval=info.err(:,1);
f2=fval;

%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,acc,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
% modify A and b
[V,D]=eig(full(A));
A_pdal=V*spdiags(sqrt(spdiags(D)),0,d,d)*V';
b_pdal=A_pdal\b;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


fmin=0.5*trueSoln'*A*trueSoln-b'*trueSoln;
%fmin=min([cy1;fval;fx;f_pdal]);
figure(2),clf,
% semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
% semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
% semilogy(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
% semilogy(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
loglog(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
loglog(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
loglog(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
loglog(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
l1=sprintf('CBCD1,%.4f s, #%d',t1,size(cy1,1)-1);
l4=sprintf('L-BFGS-B,%.4f s, #%d',t4,size(fval,1)-1);
l5=sprintf('APGD,%.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL,%.4f s, #%d',t6,size(f_pdal,1)-1);

legend(l6,l1,l4,l5);
title('Convergence of CBCD1, L-BFGS-B, APGD and PDAL');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
axis([0 9000 1E-10 3E5]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'fontsize',16,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_reorder_rcm.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_reorder_rcm.mat',...
    'f_pdal','cy1','fval','fx','fmin','t1','t4','t5','t6');
%% after amd
A=A2;b=b2;
%% CBCD
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,acc,0,1,init);
t1=toc(t0);

%% LBFGS-B
addpath /home/yu/LBSGF
addpath /home/yu/LBSGF/L-BFGS-B-C/
addpath /home/yu/LBSGF/L-BFGS-B-C/Matlab/
addpath /home/yu/LBSGF/L-BFGS-B-C/src/
disp('=== My test function, 2D === ');
n = d;
f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;
% constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;
% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
%trueSoln = xmin(p2);
trueSoln = cx1;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = acc;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 20000;
% The {f,g} is another way to call it
t0=tic;
[x,f,info] = lbfgsb( {f,g} , l, u, opts );
t4=toc(t0);

if abs(f) < 1e-8
    disp('Success!');
% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
else
    disp('Something didn''t work right :-(  ');
end
fval=info.err(:,1);
f3=fval;

%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,acc,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
% modify A and b
[V,D]=eig(full(A));
A_pdal=V*spdiags(sqrt(spdiags(D)),0,d,d)*V';
b_pdal=A_pdal\b;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters, acc,init);
t6=toc(t0);


fmin=0.5*trueSoln'*A*trueSoln-b'*trueSoln;
%fmin=min([cy1;fval;fx;f_pdal]);
figure(3),clf,
% semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
% semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
% semilogy(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
% semilogy(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
loglog(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
loglog(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
loglog(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
loglog(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
l1=sprintf('CBCD1,%.4f s, #%d',t1,size(cy1,1)-1);
l4=sprintf('L-BFGS-B,%.4f s, #%d',t4,size(fval,1)-1);
l5=sprintf('APGD,%.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL,%.4f s, #%d',t6,size(f_pdal,1)-1);

legend(l6,l1,l4,l5);
title('Convergence of CBCD1, L-BFGS-B, APGD and PDAL');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
axis([0 9000 1E-10 3E5]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'fontsize',16,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_reorder_amd.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_reorder_amd.mat',...
    'f_pdal','cy1','fval','fx','fmin','t1','t4','t5','t6');




%% type I case
% A & b
rng(1);

load /home/yu/datasets/SNAP/Social_networks/mat/facebook_combined.mat; iters=2000;
%load /home/yu/datasets/SNAP/Collaboration_networks/mat/ca-CondMat.mat; iters=400;
%load /home/yu/datasets/SNAP/Social_networks/mat/gplus_combined.mat; iters=2000;
d=size(A,1);
% d=4000;
% e1=ones(d,1);
% A=spdiags([-e1,2*e1,-e1],[-1,0,1],d,d);
% p=randperm(d,1);A=A(p,p);
A=A+speye(d)*1E-4;
%xmin=rand(d,1)*2-0.5;
xmin=rand(d,1)*0.5+0.25;
b=A*xmin;
init=0.5;
acc=1E-10;
% reorder
p1=symamd(A);
A1=A(p1,p1);b1=b(p1);
p2=symrcm(A);
A2=A(p2,p2);b2=b(p2);
% solve 


f_pdal=0;
t6=0;

%% CBCD
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,acc,0,1,init);
t1=toc(t0);

%% LBFGS-B
addpath /home/yu/LBSGF
addpath /home/yu/LBSGF/L-BFGS-B-C/
addpath /home/yu/LBSGF/L-BFGS-B-C/Matlab/
addpath /home/yu/LBSGF/L-BFGS-B-C/src/
disp('=== My test function, 2D === ');
n = d;
f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;
% constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;
% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin;
%trueSoln = cx1;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = acc;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 20000;
% The {f,g} is another way to call it
t0=tic;
[x,f,info] = lbfgsb( {f,g} , l, u, opts );
t4=toc(t0);

if abs(f) < 1e-8
    disp('Success!');
% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
else
    disp('Something didn''t work right :-(  ');
end
fval=info.err(:,1);
f1=fval;

%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,acc,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
% modify A and b
[V,D]=eig(full(A));
A_pdal=V*spdiags(sqrt(spdiags(D)),0,d,d)*V';
b_pdal=A_pdal\b;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


fmin=0.5*trueSoln'*A*trueSoln-b'*trueSoln;
%fmin=min([cy1;fval;fx;f_pdal]);
figure(4),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
semilogy(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
semilogy(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
% loglog(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
% loglog(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
% loglog(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
% loglog(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
l1=sprintf('CBCD1,%.4f s, #%d',t1,size(cy1,1)-1);
l4=sprintf('L-BFGS-B,%.4f s, #%d',t4,size(fval,1)-1);
l5=sprintf('APGD,%.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL,%.4f s, #%d',t6,size(f_pdal,1)-1);

legend(l6,l1,l4,l5);
title('Convergence of CBCD1, L-BFGS-B, APGD and PDAL');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
axis([0 2000 1E-10 1E4]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'fontsize',16,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_1reorder0.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_1reorder0.mat',...
    'f_pdal','cy1','fval','fx','fmin','t1','t4','t5','t6');

%% after rcm
A=A1;b=b1;
%% CBCD
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,acc,0,1,init);
t1=toc(t0);

%% LBFGS-B
addpath /home/yu/LBSGF
addpath /home/yu/LBSGF/L-BFGS-B-C/
addpath /home/yu/LBSGF/L-BFGS-B-C/Matlab/
addpath /home/yu/LBSGF/L-BFGS-B-C/src/
disp('=== My test function, 2D === ');
n = d;
f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;
% constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;
% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin(p1);
%trueSoln = cx1;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = acc;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 20000;
% The {f,g} is another way to call it
t0=tic;
[x,f,info] = lbfgsb( {f,g} , l, u, opts );
t4=toc(t0);

if abs(f) < 1e-8
    disp('Success!');
% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
else
    disp('Something didn''t work right :-(  ');
end
fval=info.err(:,1);
f2=fval;

%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,acc,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
% modify A and b
[V,D]=eig(full(A));
A_pdal=V*spdiags(sqrt(spdiags(D)),0,d,d)*V';
b_pdal=A_pdal\b;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


fmin=0.5*trueSoln'*A*trueSoln-b'*trueSoln;
%fmin=min([cy1;fval;fx;f_pdal]);
figure(5),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
semilogy(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
semilogy(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
% loglog(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
% loglog(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
% loglog(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
% loglog(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
l1=sprintf('CBCD1,%.4f s, #%d',t1,size(cy1,1)-1);
l4=sprintf('L-BFGS-B,%.4f s, #%d',t4,size(fval,1)-1);
l5=sprintf('APGD,%.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL,%.4f s, #%d',t6,size(f_pdal,1)-1);

legend(l6,l1,l4,l5);
title('Convergence of CBCD1, L-BFGS-B, APGD and PDAL');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
axis([0 2000 1E-10 1E4]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'fontsize',16,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_1reorder_rcm.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_1reorder_rcm.mat',...
    'f_pdal','cy1','fval','fx','fmin','t1','t4','t5','t6');

%% after amd
A=A2;b=b2;
%% CBCD
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,acc,0,1,init);
t1=toc(t0);

%% LBFGS-B
addpath /home/yu/LBSGF
addpath /home/yu/LBSGF/L-BFGS-B-C/
addpath /home/yu/LBSGF/L-BFGS-B-C/Matlab/
addpath /home/yu/LBSGF/L-BFGS-B-C/src/
disp('=== My test function, 2D === ');
n = d;
f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;
% constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;
% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin(p2);
%trueSoln = cx1;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = acc;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 20000;
% The {f,g} is another way to call it
t0=tic;
[x,f,info] = lbfgsb( {f,g} , l, u, opts );
t4=toc(t0);

if abs(f) < 1e-8
    disp('Success!');
% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
else
    disp('Something didn''t work right :-(  ');
end
fval=info.err(:,1);
f3=fval;

%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,acc,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
% modify A and b
[V,D]=eig(full(A));
A_pdal=V*spdiags(sqrt(spdiags(D)),0,d,d)*V';
b_pdal=A_pdal\b;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters, acc,init);
t6=toc(t0);


fmin=0.5*trueSoln'*A*trueSoln-b'*trueSoln;
%fmin=min([cy1;fval;fx;f_pdal]);
figure(6),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
semilogy(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
semilogy(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
% loglog(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',2);hold on;
% loglog(1:size(cy1,1),cy1-fmin,'r','LineWidth',2);hold on;
% loglog(1:size(fval,1),fval-fmin,'c','LineWidth',2);hold on;
% loglog(1:size(fx,1),fx-fmin,'m','LineWidth',2);hold on;
l1=sprintf('CBCD1,%.4f s, #%d',t1,size(cy1,1)-1);
l4=sprintf('L-BFGS-B,%.4f s, #%d',t4,size(fval,1)-1);
l5=sprintf('APGD,%.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL,%.4f s, #%d',t6,size(f_pdal,1)-1);

legend(l6,l1,l4,l5);
title('Convergence of CBCD1, L-BFGS-B, APGD and PDAL');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
axis([0 2000 1E-10 1E4]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'fontsize',16,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_1reorder_amd.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_reorder_4methods/chapter_1reorder_amd.mat',...
    'f_pdal','cy1','fval','fx','fmin','t1','t4','t5','t6');






