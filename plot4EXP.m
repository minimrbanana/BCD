function plot4EXP(idx,EXP)
%% plot from EXP.mat
% if EXP is given, as called in the main_exp function, use it, 
% else as called alone, load EXP from the 'result' folder
close all;
% dir to save the plots
saveDir = ['./img/EXP_idx_' num2str(idx) '/'];
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

if nargin ==1
    dir = ['./result/EXP_idx_' num2str(idx) '/EXP.mat'];
    load(dir);
end
%% plot
% see if plot the matrix
if EXP.isplot == 1
    figure(1),clf;
    h1 = imagesc(EXP.A);
    colorbar; hold all;
    title('Matrix A');
    saveas(h1,[EXP.output_dir 'figure1.png']);
    figure(2),clf;
    h2 = imagesc(EXP.B);
    colorbar; hold all;
    title('Matrix B');
    saveas(h2,[EXP.output_dir 'figure2.png']);
    figure(3),clf;
    h3 = imagesc(EXP.C);
    colorbar; hold all;
    title('Matrix C');
    saveas(h3,[EXP.output_dir 'figure3.png']);
end
%% if plot the convergence
% first choose the index of experiment with 
%               %
% mean of iters %
%               %
% number of epochs for each experiment with mat A and CBCD1
% is read from EXP.epoch1
% sort the number of epochs and choose the middle one
[~,I] = sort(EXP.epoch1(1:EXP.n_loop));
index = I(ceil(EXP.n_loop/2));
%%
% then plot the cooresponding curve of convergence
figure(4),clf;
% in figure 4 we show the convergence of the KKT condition
% which is based on the normal cone

% plot matrix A
semilogy(0:size(EXP.KKT_A{index},1)-1,...
    EXP.KKT_A{index},'r','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKT_A{index+EXP.n_loop},1)-1,...
    EXP.KKT_A{index+EXP.n_loop},'g','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKT_A{index+EXP.n_loop*2},1)-1,...
    EXP.KKT_A{index+EXP.n_loop*2},'b','LineWidth',2.5);hold on;
% plot matrix A of RBCD
semilogy(0:size(EXP.KKTr_A{index},1)-1,...
    EXP.KKTr_A{index},'r--','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKTr_A{index+EXP.n_loop},1)-1,...
    EXP.KKTr_A{index+EXP.n_loop},'g--','LineWidth',2);hold on;
semilogy(0:size(EXP.KKTr_A{index+EXP.n_loop*2},1)-1,...
    EXP.KKTr_A{index+EXP.n_loop*2},'b--','LineWidth',1.5);hold on;
% legend of A: CBCD & RBCD
% first add runtime into legend
l1=sprintf('CBCD1, %.4f s, #%d',mean(EXP.T_c(1,:)),...
    size(EXP.KKT_A{index},1)-1);
l2=sprintf('CBCD2, %.4f s, #%d',mean(EXP.T_c(2,:)),...
    size(EXP.KKT_A{index+EXP.n_loop},1)-1);
l3=sprintf('CBCD3, %.4f s, #%d',mean(EXP.T_c(3,:)),...
    size(EXP.KKT_A{index+EXP.n_loop*2},1)-1);
l1r=sprintf('RBCD1, %.4f s, #%d',mean(EXP.T_r(1,:)),...
    size(EXP.KKTr_A{index},1)-1);
l2r=sprintf('RBCD2, %.4f s, #%d',mean(EXP.T_r(2,:)),...
    size(EXP.KKTr_A{index+EXP.n_loop},1)-1);
l3r=sprintf('RBCD3, %.4f s, #%d',mean(EXP.T_r(3,:)),...
    size(EXP.KKTr_A{index+EXP.n_loop*2},1)-1);
legend(l1,l2,l3,l1r,l2r,l3r);
set(gca,'fontsize',14);
xlabel('#epoch');ylabel('KKT Condition');
title(['Convergence Speed A; #' num2str(idx) '; #' num2str(index) ]);
saveas(gca,[saveDir 'figure4.png']);

%%
figure(5),clf;
% in figure 5 we show the convergence of the KKT condition
% which is based on the normal cone

% plot matrix B
semilogy(0:size(EXP.KKT_B{index},1)-1,...
    EXP.KKT_B{index},'r','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKT_B{index+EXP.n_loop},1)-1,...
    EXP.KKT_B{index+EXP.n_loop},'g','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKT_B{index+EXP.n_loop*2},1)-1,...
    EXP.KKT_B{index+EXP.n_loop*2},'b','LineWidth',2.5);hold on;
% plot matrix B of RBCD
semilogy(0:size(EXP.KKTr_B{index},1)-1,...
    EXP.KKTr_B{index},'r--','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKTr_B{index+EXP.n_loop},1)-1,...
    EXP.KKTr_B{index+EXP.n_loop},'g--','LineWidth',2);hold on;
semilogy(0:size(EXP.KKTr_B{index+EXP.n_loop*2},1)-1,...
    EXP.KKTr_B{index+EXP.n_loop*2},'b--','LineWidth',1.5);hold on;
% legend of B; CBCD & RBCD
% first add runtime into legend
l4 =sprintf('CBCD1p, %.4f s, #%d',mean(EXP.T_c(4,:)),...
    size(EXP.KKT_B{index},1)-1);
l5 =sprintf('CBCD2p, %.4f s, #%d',mean(EXP.T_c(5,:)),...
    size(EXP.KKT_B{index+EXP.n_loop},1)-1);
l6 =sprintf('CBCD3p, %.4f s, #%d',mean(EXP.T_c(6,:)),...
    size(EXP.KKT_B{index+EXP.n_loop*2},1)-1);
l4r=sprintf('RBCD1p, %.4f s, #%d',mean(EXP.T_r(4,:)),...
    size(EXP.KKTr_B{index},1)-1);
l5r=sprintf('RBCD2p, %.4f s, #%d',mean(EXP.T_r(5,:)),...
    size(EXP.KKTr_B{index+EXP.n_loop},1)-1);
l6r=sprintf('RBCD3p, %.4f s, #%d',mean(EXP.T_r(6,:)),...
    size(EXP.KKTr_B{index+EXP.n_loop*2},1)-1);
legend(l4,l5,l6,l4r,l5r,l6r);
set(gca,'fontsize',14);
xlabel('#epoch');ylabel('KKT Condition');
title(['Convergence Speed B; #EXP ' num2str(idx) '; #' num2str(index) ]);
saveas(gca,[saveDir 'figure5.png']);
%%
figure(6),clf;
% in figure 6 we show the convergence of the KKT condition
% which is based on the normal cone

% plot matrix C, first cyclic and then random
semilogy(0:size(EXP.KKT_C{index},1)-1,...
    EXP.KKT_C{index},'r','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKT_C{index+EXP.n_loop},1)-1,...
    EXP.KKT_C{index+EXP.n_loop},'g','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKT_C{index+EXP.n_loop*2},1)-1,...
    EXP.KKT_C{index+EXP.n_loop*2},'b','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKTr_C{index},1)-1,...
    EXP.KKTr_C{index},'r--','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKTr_C{index+EXP.n_loop},1)-1,...
    EXP.KKTr_C{index+EXP.n_loop},'g--','LineWidth',2.5);hold on;
semilogy(0:size(EXP.KKTr_C{index+EXP.n_loop*2},1)-1,...
    EXP.KKTr_C{index+EXP.n_loop*2},'b--','LineWidth',2.5);hold on;
% legend of C
% first add runtime into legend
l7 =sprintf('CBCD1r,   %.4f s, #%d',mean(EXP.T_c(7,:)),...
    size(EXP.KKT_C{index},1)-1);
l8 =sprintf('CBCD2r,   %.4f s, #%d',mean(EXP.T_c(8,:)),...
    size(EXP.KKT_C{index+EXP.n_loop},1)-1);
l9 =sprintf('CBCD3r,   %.4f s, #%d',mean(EXP.T_c(9,:)),...
    size(EXP.KKT_C{index+EXP.n_loop*2},1)-1);
l7r=sprintf('RBCD1r,   %.4f s, #%d',mean(EXP.T_r(7,:)),...
    size(EXP.KKTr_C{index},1)-1);
l8r=sprintf('RBCD2r,   %.4f s, #%d',mean(EXP.T_r(8,:)),...
    size(EXP.KKTr_C{index+EXP.n_loop},1)-1);
l9r=sprintf('RBCD3r,   %.4f s, #%d',mean(EXP.T_r(9,:)),...
    size(EXP.KKTr_C{index+EXP.n_loop*2},1)-1);
legend(l7,l8,l9,l7r,l8r,l9r);
set(gca,'fontsize',14);
xlabel('#epoch');ylabel('KKT Condition');
RCMstr = sprintf(';RCM=%.4fs',EXP.tRCM);
title(['Convergence Speed C; #' num2str(idx) '-#' num2str(index) RCMstr]);
saveas(gca,[saveDir 'figure6.png']);

%%
% first get all the function values according to the index
% also get the min of f(x) from these values
fx1a = EXP.OBJ_A{index};%fx1a = fx1a(1:nnz(fx1a)+1);
fx2a = EXP.OBJ_A{index+EXP.n_loop};%fx2a = fx2a(1:nnz(fx2a)+1);
fx3a = EXP.OBJ_A{index+EXP.n_loop*2};%fx3a = fx3a(1:nnz(fx3a)+1);
fx1b = EXP.OBJ_B{index};%fx1b = fx1b(1:nnz(fx1b)+1);
fx2b = EXP.OBJ_B{index+EXP.n_loop};%fx2b = fx2b(1:nnz(fx2b)+1);
fx3b = EXP.OBJ_B{index+EXP.n_loop*2};%fx3b = fx3b(1:nnz(fx3b)+1);
fx1c = EXP.OBJ_C{index};%fx1c = fx1c(1:nnz(fx1c)+1);
fx2c = EXP.OBJ_C{index+EXP.n_loop};%fx2c = fx2c(1:nnz(fx2c)+1);
fx3c = EXP.OBJ_C{index+EXP.n_loop*2};%fx3c = fx3c(1:nnz(fx3c)+1);
fmin = min([fx1a;fx2a;fx3a;fx1b;fx2b;fx3b;fx1c;fx2c;fx3c]);
%% then plot function value for matrix A
figure(7),clf;
% in figure 7 we show the convergence of the Objective
% which is the function value f(x^k)
% plot matrix A
semilogy(0:size(fx1a,1)-1,fx1a-fmin+eps,...
    'r','LineWidth',2.5);hold on;
semilogy(0:size(fx2a,1)-1,fx2a-fmin+eps,...
    'g','LineWidth',2.5);hold on;
semilogy(0:size(fx3a,1)-1,fx3a-fmin+eps,...
    'b','LineWidth',2.5);hold on;
% legend of A & B
l1=sprintf('CBCD1,    #%d',size(fx1a,1)-1);
l2=sprintf('CBCD2,    #%d',size(fx2a,1)-1);
l3=sprintf('CBCD3,    #%d',size(fx3a,1)-1);
legend(l1,l2,l3);
set(gca,'fontsize',14);
xlabel('#epoch');ylabel('f(x^k)-p^*+eps');
title(['Convergence Speed mat A #EXP ' num2str(idx) ' #' num2str(index) ]);
saveas(gca,[saveDir 'figure7.png']);

%% then plot function value for matrix B
figure(8),clf;
% in figure 8 we show the convergence of the Objective
% which is the function value f(x^k)
% plot matrix B
semilogy(0:size(fx1b,1)-1,fx1b-fmin+eps,...
    'r','LineWidth',2.5);hold on;
semilogy(0:size(fx2b,1)-1,fx2b-fmin+eps,...
    'g','LineWidth',2.5);hold on;
semilogy(0:size(fx3b,1)-1,fx3b-fmin+eps,...
    'b','LineWidth',2.5);hold on;
% legend of A & B
l4=sprintf('CBCD1p,  #%d',size(fx1b,1)-1);
l5=sprintf('CBCD2p,  #%d',size(fx2b,1)-1);
l6=sprintf('CBCD3p,  #%d',size(fx3b,1)-1);
legend(l4,l5,l6);
set(gca,'fontsize',14);
xlabel('#epoch');ylabel('f(x^k)-p^*+eps');
title(['Convergence Speed mat B #EXP ' num2str(idx) ' #' num2str(index) ]);
saveas(gca,[saveDir 'figure8.png']);

%%
figure(9),clf;
% in figure 9 we show the convergence of the Objective
% which is the function value f(x^k)
% plot matrix C
semilogy(0:size(fx1c,1)-1,fx1c-fmin+eps,'r','LineWidth',2.5);hold on;
semilogy(0:size(fx2c,1)-1,fx2c-fmin+eps,'g','LineWidth',2.5);hold on;
semilogy(0:size(fx3c,1)-1,fx3c-fmin+eps,'b','LineWidth',2.5);hold on;
% legend of A & B
l7=sprintf('CBCD1r,    #%d',size(fx1c,1)-1);
l8=sprintf('CBCD2r,    #%d',size(fx2c,1)-1);
l9=sprintf('CBCD3r,    #%d',size(fx3c,1)-1);
legend(l7,l8,l9);
set(gca,'fontsize',14);
xlabel('#epoch');ylabel('f(x^k)-p^*+eps');
title(['Convergence Speed mat C #EXP ' num2str(idx) ' #' num2str(index) ]);
saveas(gca,[saveDir 'figure9.png']);


end
