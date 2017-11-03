%% compare the 4 methods on facebook matrix, which has band structure
addpath /home/yu/bcd/BCD/
addpath /home/yu/LBSGF
addpath /home/yu/LBSGF/L-BFGS-B-C/
addpath /home/yu/LBSGF/L-BFGS-B-C/Matlab/
addpath /home/yu/LBSGF/L-BFGS-B-C/src/
load /home/yu/datasets/SNAP/Social_networks/mat/facebook_combined.mat
%%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% HAVE TO MODIFY F_MIN AND X_MIN WHEN CHANGING CONSTRAINTS
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%
d=size(A,1);
diag=sum(A);
A=A+speye(d)*0.5;
rng(1);
%% Type I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    type I , init = 0.5     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin=rand(d,1)*0.5+0.25;
b=A*xmin;

%% true minimum and f min
fmin=0.5*xmin'*A*xmin-b'*xmin;
init=0.5;
Fsize=17;
%% CBCD solver
iters=2000;
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,1E-10,0,1,init);
t1=toc(t0);
t0=tic;
[cx2, cy2] = CBCD_size2_fx(A, b, d, iters,1E-10,0,1,init);
t2=toc(t0);
t0=tic;
[cx3, cy3] = CBCD_size3_fx(A, b, d, iters,1E-10,0,1,init);
t3=toc(t0);

%% lbfgs-b solver
disp('=== My test function, 2D === ');
n = d;

f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;

% There are no constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 10000;
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
f_bfgs=info.err(:,1);


%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,1E-10,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
% modify A and b
[V,D]=eig(full(A));
A_pdal=V*spdiags(sqrt(spdiags(D)),0,d,d)*V';
b_pdal=A_pdal\b;
acc=1E-10;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


%% plot
figure(1),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',4);
hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',4);
hold on;
semilogy(1:size(cy2,1),cy2-fmin,'g','LineWidth',3);
hold on;
semilogy(1:size(cy3,1),cy3-fmin,'b','LineWidth',2);
hold on;
semilogy(1:length(f_bfgs),f_bfgs-fmin,'c','LineWidth',2);
hold on;
semilogy(1:length(fx),fx-fmin,'m','LineWidth',2);
hold on;


l1=sprintf('CBCD1,   %.4f s, #%d',t1,size(cy1,1)-1);
l2=sprintf('CBCD2,   %.4f s, #%d',t2,size(cy2,1)-1);
l3=sprintf('CBCD3,   %.4f s, #%d',t3,size(cy3,1)-1);
l4=sprintf('L-BFGS-B,%.2f s, #%d',t4,size(f_bfgs,1)-1);
l5=sprintf('AGD  ,   %.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL ,   %.4f s, #%d',t6,size(f_pdal,1)'-1);

legend(l6,l1,l2,l3,l4,l5);
title('Convergence of CBCD, L-BFGS-B, PDAL and APGD');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_5.png');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_5.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_5.mat',...
    'f_pdal','cy1','cy2','cy3','f_bfgs','fx','fmin','t1','t2','t3','t4','t5','t6');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    type I , init = 0.0     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% true minimum and f min
init=0;

%% CBCD solver
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,1E-10,0,1,init);
t1=toc(t0);
t0=tic;
[cx2, cy2] = CBCD_size2_fx(A, b, d, iters,1E-10,0,1,init);
t2=toc(t0);
t0=tic;
[cx3, cy3] = CBCD_size3_fx(A, b, d, iters,1E-10,0,1,init);
t3=toc(t0);

%% lbfgs-b solver
disp('=== My test function, 2D === ');
n = d;

f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;

% There are no constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 10000;
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
f_bfgs=info.err(:,1);


%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,1E-10,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


%% plot
figure(2),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',4);
hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',4);
hold on;
semilogy(1:size(cy2,1),cy2-fmin,'g','LineWidth',3);
hold on;
semilogy(1:size(cy3,1),cy3-fmin,'b','LineWidth',2);
hold on;
semilogy(1:length(f_bfgs),f_bfgs-fmin,'c','LineWidth',2);
hold on;
semilogy(1:length(fx),fx-fmin,'m','LineWidth',2);
hold on;


l1=sprintf('CBCD1,   %.4f s, #%d',t1,size(cy1,1)-1);
l2=sprintf('CBCD2,   %.4f s, #%d',t2,size(cy2,1)-1);
l3=sprintf('CBCD3,   %.4f s, #%d',t3,size(cy3,1)-1);
l4=sprintf('L-BFGS-B,%.2f s, #%d',t4,size(f_bfgs,1)-1);
l5=sprintf('AGD  ,   %.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL ,   %.4f s, #%d',t6,size(f_pdal,1)'-1);

legend(l6,l1,l2,l3,l4,l5);
title('Convergence of CBCD, L-BFGS-B, PDAL and APGD');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_0.png');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_0.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_0.mat',...
    'f_pdal','cy1','cy2','cy3','f_bfgs','fx','fmin','t1','t2','t3','t4','t5','t6');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    type I , init = 1.0     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% true minimum and f min
init=1;

%% CBCD solver
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,1E-10,0,1,init);
t1=toc(t0);
t0=tic;
[cx2, cy2] = CBCD_size2_fx(A, b, d, iters,1E-10,0,1,init);
t2=toc(t0);
t0=tic;
[cx3, cy3] = CBCD_size3_fx(A, b, d, iters,1E-10,0,1,init);
t3=toc(t0);

%% lbfgs-b solver
disp('=== My test function, 2D === ');
n = d;

f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;

% There are no constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 10000;
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
f_bfgs=info.err(:,1);


%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,1E-10,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


%% plot
figure(3),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',4);
hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',4);
hold on;
semilogy(1:size(cy2,1),cy2-fmin,'g','LineWidth',3);
hold on;
semilogy(1:size(cy3,1),cy3-fmin,'b','LineWidth',2);
hold on;
semilogy(1:length(f_bfgs),f_bfgs-fmin,'c','LineWidth',2);
hold on;
semilogy(1:length(fx),fx-fmin,'m','LineWidth',2);
hold on;


l1=sprintf('CBCD1,   %.4f s, #%d',t1,size(cy1,1)-1);
l2=sprintf('CBCD2,   %.4f s, #%d',t2,size(cy2,1)-1);
l3=sprintf('CBCD3,   %.4f s, #%d',t3,size(cy3,1)-1);
l4=sprintf('L-BFGS-B,%.2f s, #%d',t4,size(f_bfgs,1)-1);
l5=sprintf('AGD  ,   %.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL ,   %.4f s, #%d',t6,size(f_pdal,1)'-1);

legend(l6,l1,l2,l3,l4,l5);
title('Convergence of CBCD, L-BFGS-B, PDAL and APGD');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_1.png');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_1.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/uNest_1.mat',...
    'f_pdal','cy1','cy2','cy3','f_bfgs','fx','fmin','t1','t2','t3','t4','t5','t6');
%% Type II  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    type II , init = 0.5     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin0=rand(d,1)*2-0.5;
b=A*xmin0;

%% true minimum and f min
init=0.5;

%% CBCD solver
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,1E-10,0,1,init);
t1=toc(t0);
t0=tic;
[cx2, cy2] = CBCD_size2_fx(A, b, d, iters,1E-10,0,1,init);
t2=toc(t0);
t0=tic;
[cx3, cy3] = CBCD_size3_fx(A, b, d, iters,1E-10,0,1,init);
t3=toc(t0);

%% true minimum and f min
xmin = cx1;
fmin=0.5*xmin'*A*xmin-b'*xmin;
%% lbfgs-b solver
disp('=== My test function, 2D === ');
n = d;

f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;

% There are no constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 10000;
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
f_bfgs=info.err(:,1);


%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,1E-10,iters);
%fx=[0;0;0];
t5=toc(t0);


%% Chambolle Pock primal dual
% change b here
b_pdal=A_pdal\b;
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


fmin=min([f_pdal',cy1',cy2',cy3',f_bfgs',fx']);
%% plot
figure(4),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',4);
hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',4);
hold on;
semilogy(1:size(cy2,1),cy2-fmin,'g','LineWidth',3);
hold on;
semilogy(1:size(cy3,1),cy3-fmin,'b','LineWidth',2);
hold on;
semilogy(1:length(f_bfgs),f_bfgs-fmin,'c','LineWidth',2);
hold on;
semilogy(1:length(fx),fx-fmin,'m','LineWidth',2);
hold on;


l1=sprintf('CBCD1,   %.4f s, #%d',t1,size(cy1,1)-1);
l2=sprintf('CBCD2,   %.4f s, #%d',t2,size(cy2,1)-1);
l3=sprintf('CBCD3,   %.4f s, #%d',t3,size(cy3,1)-1);
l4=sprintf('L-BFGS-B,%.2f s, #%d',t4,size(f_bfgs,1)-1);
l5=sprintf('AGD  ,   %.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL ,   %.4f s, #%d',t6,size(f_pdal,1)'-1);

legend(l6,l1,l2,l3,l4,l5);
title('Convergence of CBCD, L-BFGS-B, PDAL and APGD');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_5.png');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_5.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_5.mat',...
    'f_pdal','cy1','cy2','cy3','f_bfgs','fx','fmin','t1','t2','t3','t4','t5','t6');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    type II , init = 0.0     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% true minimum and f min
init=0;

%% CBCD solver
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,1E-10,0,1,init);
t1=toc(t0);
t0=tic;
[cx2, cy2] = CBCD_size2_fx(A, b, d, iters,1E-10,0,1,init);
t2=toc(t0);
t0=tic;
[cx3, cy3] = CBCD_size3_fx(A, b, d, iters,1E-10,0,1,init);
t3=toc(t0);

%% lbfgs-b solver
disp('=== My test function, 2D === ');
n = d;

f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;

% There are no constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 10000;
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
f_bfgs=info.err(:,1);


%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,1E-10,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


fmin=min([f_pdal',cy1',cy2',cy3',f_bfgs',fx']);
%% plot
figure(5),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',4);
hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',4);
hold on;
semilogy(1:size(cy2,1),cy2-fmin,'g','LineWidth',3);
hold on;
semilogy(1:size(cy3,1),cy3-fmin,'b','LineWidth',2);
hold on;
semilogy(1:length(f_bfgs),f_bfgs-fmin,'c','LineWidth',2);
hold on;
semilogy(1:length(fx),fx-fmin,'m','LineWidth',2);
hold on;


l1=sprintf('CBCD1,   %.4f s, #%d',t1,size(cy1,1)-1);
l2=sprintf('CBCD2,   %.4f s, #%d',t2,size(cy2,1)-1);
l3=sprintf('CBCD3,   %.4f s, #%d',t3,size(cy3,1)-1);
l4=sprintf('L-BFGS-B,%.2f s, #%d',t4,size(f_bfgs,1)-1);
l5=sprintf('AGD  ,   %.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL ,   %.4f s, #%d',t6,size(f_pdal,1)'-1);

legend(l6,l1,l2,l3,l4,l5);
title('Convergence of CBCD, L-BFGS-B, PDAL and APGD');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_0.png');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_0.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_0.mat',...
    'f_pdal','cy1','cy2','cy3','f_bfgs','fx','fmin','t1','t2','t3','t4','t5','t6');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    type I , init = 1.0     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% true minimum and f min
init=1;

%% CBCD solver
t0=tic;
[cx1, cy1] = CBCD_size1_fx(A, b, d, iters,1E-10,0,1,init);
t1=toc(t0);
t0=tic;
[cx2, cy2] = CBCD_size2_fx(A, b, d, iters,1E-10,0,1,init);
t2=toc(t0);
t0=tic;
[cx3, cy3] = CBCD_size3_fx(A, b, d, iters,1E-10,0,1,init);
t3=toc(t0);

%% lbfgs-b solver
disp('=== My test function, 2D === ');
n = d;

f = @(x) 0.5*x'*A*x-b'*x;
g = @(x) A*x-b;

% There are no constraints
l   = zeros(n,1);%-inf(n,1);
u   = ones(n,1);%inf(n,1);

opts    = struct( 'x0', ones(d,1)*init );
opts.printEvery     = 1;
opts.m  = 30;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = xmin;
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1;
opts.maxIts     = iters;
opts.maxTotalIts= 10000;
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
f_bfgs=info.err(:,1);


%% nesterov's accelerated projected gradient
addpath /home/yu/FPG
t0=tic;
[x,fx]=FPG(A,b,init,1E-10,iters);
%fx=[0;0;0];
t5=toc(t0);

%% Chambolle Pock primal dual
t0=tic;
[x_pdal, y_pdal, epoch_pdal, f_pdal] = PDAL(A_pdal, b_pdal, iters,acc, init);
t6=toc(t0);


fmin=min([f_pdal',cy1',cy2',cy3',f_bfgs',fx']);
%% plot
figure(6),clf,
semilogy(1:size(f_pdal,1),f_pdal-fmin,'color',[1,0.5,0],'LineWidth',4);
hold on;
semilogy(1:size(cy1,1),cy1-fmin,'r','LineWidth',4);
hold on;
semilogy(1:size(cy2,1),cy2-fmin,'g','LineWidth',3);
hold on;
semilogy(1:size(cy3,1),cy3-fmin,'b','LineWidth',2);
hold on;
semilogy(1:length(f_bfgs),f_bfgs-fmin,'c','LineWidth',2);
hold on;
semilogy(1:length(fx),fx-fmin,'m','LineWidth',2);
hold on;


l1=sprintf('CBCD1,   %.4f s, #%d',t1,size(cy1,1)-1);
l2=sprintf('CBCD2,   %.4f s, #%d',t2,size(cy2,1)-1);
l3=sprintf('CBCD3,   %.4f s, #%d',t3,size(cy3,1)-1);
l4=sprintf('L-BFGS-B,%.2f s, #%d',t4,size(f_bfgs,1)-1);
l5=sprintf('AGD  ,   %.4f s, #%d',t5,size(fx,1)-1);
l6=sprintf('PDAL ,   %.4f s, #%d',t6,size(f_pdal,1)'-1);

legend(l6,l1,l2,l3,l4,l5);
title('Convergence of CBCD, L-BFGS-B, PDAL and APGD');
xlabel('#epoch');
ylabel('log(f(x^k)-f^*)');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_1.png');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_1.pdf');
save('/home/yu/bcd/BCD/test_scripts/thesis_compare/Facebook/Nest_1.mat',...
    'f_pdal','cy1','cy2','cy3','f_bfgs','fx','fmin','t1','t2','t3','t4','t5','t6');