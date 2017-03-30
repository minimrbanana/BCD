function main_exp(eidx)
% the main function for the experiments and save the results
% clear and then load the experimental parameters
close all;
clear EXP;
EXP = exp_detail(eidx);
d = EXP.d;
iters = EXP.max_iter;
% do permutation on A, we get B, reordering on B we get C
perm = randperm(EXP.d);
B = EXP.A(:,perm);
B = B(perm,:);
%[C,reod] = RCM(B);
reod = symrcm(B);
C = B(reod,reod);
% see if plot
if EXP.isplot == 1
    figure(1),clf;
    h1 = imagesc(EXP.A);
    colorbar; hold all;
    title('Matrix A');
    saveas(h1,[EXP.output_dir 'figure1.png']);
    figure(2),clf;
    h2 = imagesc(B);
    colorbar; hold all;
    title('Matrix B');
    saveas(h2,[EXP.output_dir 'figure2.png']);
    figure(3),clf;
    h3 = imagesc(C);
    colorbar; hold all;
    title('Matrix C');
    saveas(h3,[EXP.output_dir 'figure3.png']);
end
% create variables for saving the results
% save runtime
T = zeros(9,EXP.n_loop);
% save the objective, related to normal cone
% 1,2,3 for matrix A, B, C
Obj1 = zeros(3*EXP.n_loop,iters);
% first n_loop rows store size1
% second n_loop rows store size2
% third n_loop rows store size3
%
Obj2 = zeros(3*EXP.n_loop,iters);
Obj3 = zeros(3*EXP.n_loop,iters);
% loop to average the cenvergence
for loop=1:EXP.n_loop
    b = randn(EXP.d,1);
    Bb= b(perm);
    Cb = Bb(reod);
    %
    % runtime of matrix A
    %
    profile on;
    [~, y1a] = CBCD_size1_mex_sparse(EXP.A, b, d, iters);
    [~, y2a] = CBCD_size2_ss(EXP.A, b, d, iters);
    [~, y3a] = CBCD_size3_ss(EXP.A, b, d, iters);
    p=profile('info');
    profile off;
    % save time for matrix A
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_mex_sparse')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(1,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(1,loop));
            fprintf('#epochs : %d \n',length(y1a));
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(2,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(2,loop));
            fprintf('#epochs : %d \n',length(y2a));
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(3,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(3,loop));
            fprintf('#epochs : %d \n',length(y3a));
        end
    end
    % save objective
    Obj1(loop,1:length(y1a)) = y1a;
    Obj1(loop+EXP.n_loop,1:length(y2a)) = y2a;
    Obj1(loop+EXP.n_loop*2,1:length(y3a)) = y3a;
    %
    % runtime of matrix B
    %
    profile on;
    [~, y1b] = CBCD_size1_mex_sparse(B, Bb, d, iters);
    [~, y2b] = CBCD_size2_ss(B, Bb, d, iters);
    [~, y3b] = CBCD_size3_ss(B, Bb, d, iters);
    p=profile('info');
    profile off;
    % save time for matrix B
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_mex_sparse')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(4,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(4,loop));
            fprintf('#epochs : %d \n',length(y1b));
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(5,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(5,loop));
            fprintf('#epochs : %d \n',length(y2b));
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(6,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(6,loop));
            fprintf('#epochs : %d \n',length(y3b));
        end
    end
    % save objective
    Obj2(loop,1:length(y1b)) = y1b;
    Obj2(loop+EXP.n_loop,1:length(y2b)) = y2b;
    Obj2(loop+EXP.n_loop*2,1:length(y3b)) = y3b;
    %
    % runtime of matrix C
    %
    profile on;
    [~, y1c] = CBCD_size1_mex_sparse(C, Cb, d, iters);
    [~, y2c] = CBCD_size2_ss(C, Cb, d, iters);
    [~, y3c] = CBCD_size3_ss(C, Cb, d, iters);
    p=profile('info');
    profile off;
    % save time for matrix C
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_mex_sparse')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(7,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(7,loop));
            fprintf('#epochs : %d \n',length(y1c));
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(8,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(8,loop));
            fprintf('#epochs : %d \n',length(y2c));
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(9,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(9,loop));
            fprintf('#epochs : %d \n',length(y3c));
        end
    end
    % save objective
    Obj3(loop,1:length(y1c)) = y1c;
    Obj3(loop+EXP.n_loop,1:length(y2c)) = y2c;
    Obj3(loop+EXP.n_loop*2,1:length(y3c)) = y3c;
end
EXP.T = T;
EXP.Obj1 = Obj1;
EXP.Obj2 = Obj2;
EXP.Obj3 = Obj3;
EXP.B = B;
EXP.C = C;
if EXP.save==1
    save([EXP.output_dir 'EXP.mat'],'EXP');
end

% if plot the convergence
if EXP.plot_convergence ==1
    figure(4),clf;
    % plot matrix A
    semilogy(0:EXP.max_iter-1,median(EXP.Obj1(1:EXP.n_loop,:),1),'r','LineWidth',2.5);hold on;
    semilogy(0:EXP.max_iter-1,median(EXP.Obj1(1+EXP.n_loop:EXP.n_loop*2,:),1),'g','LineWidth',2.5);hold on;
    semilogy(0:EXP.max_iter-1,median(EXP.Obj1(1+EXP.n_loop*2:end,:),1),'b','LineWidth',2.5);hold on;
    % plot matrix B
    semilogy(0:EXP.max_iter-1,median(EXP.Obj2(1:EXP.n_loop,:),1),'r--','LineWidth',2.5);hold on;
    semilogy(0:EXP.max_iter-1,median(EXP.Obj2(1+EXP.n_loop:EXP.n_loop*2,:),1),'g--','LineWidth',2);hold on;
    semilogy(0:EXP.max_iter-1,median(EXP.Obj2(1+EXP.n_loop*2:end,:),1),'b--','LineWidth',1.5);hold on;
    % legend of A & B
    % first add runtime into legend
    l1=sprintf('CBCD1,   %.4f s, #%d',mean(EXP.T(1,:)),nnz(median(EXP.Obj1(1:EXP.n_loop,:),1)));
    l2=sprintf('CBCD2,   %.4f s, #%d',mean(EXP.T(2,:)),nnz(median(EXP.Obj1(1+EXP.n_loop:EXP.n_loop*2,:),1)));
    l3=sprintf('CBCD3,   %.4f s, #%d',mean(EXP.T(3,:)),nnz(median(EXP.Obj1(1+EXP.n_loop*2:end,:),1)));
    l4=sprintf('CBCD1p, %.4f s, #%d',mean(EXP.T(4,:)),nnz(median(EXP.Obj2(1:EXP.n_loop,:),1)));
    l5=sprintf('CBCD2p, %.4f s, #%d',mean(EXP.T(5,:)),nnz(median(EXP.Obj2(1+EXP.n_loop:EXP.n_loop*2,:),1)));
    l6=sprintf('CBCD3p, %.4f s, #%d',mean(EXP.T(6,:)),nnz(median(EXP.Obj2(1+EXP.n_loop*2:end,:),1)));
    legend(l1,l2,l3,l4,l5,l6);
    xlabel('#epoch');ylabel('Objective');
    title(['Convergence Speed and Runtime mat A&B #EXP ' num2str(eidx)]);
    set(gca,'fontsize',14);
    saveas(gca,[EXP.output_dir 'figure4.png']);

    figure(5),clf;
    % plot matrix C
    semilogy(0:EXP.max_iter-1,median(EXP.Obj3(1:EXP.n_loop,:),1),'r','LineWidth',2.5);hold on;
    semilogy(0:EXP.max_iter-1,median(EXP.Obj3(1+EXP.n_loop:EXP.n_loop*2,:),1),'g','LineWidth',2.5);hold on;
    semilogy(0:EXP.max_iter-1,median(EXP.Obj3(1+EXP.n_loop*2:end,:),1),'b','LineWidth',2.5);hold on;
    % legend of C
    % first add runtime into legend
    l7=sprintf('CBCD1r,   %.4f s, #%d',mean(EXP.T(7,:)),nnz(median(EXP.Obj3(1:EXP.n_loop,:),1)));
    l8=sprintf('CBCD2r,   %.4f s, #%d',mean(EXP.T(8,:)),nnz(median(EXP.Obj3(1+EXP.n_loop:EXP.n_loop*2,:),1)));
    l9=sprintf('CBCD3r,   %.4f s, #%d',mean(EXP.T(9,:)),nnz(median(EXP.Obj3(1+EXP.n_loop*2:end,:),1)));
    legend(l7,l8,l9);
    xlabel('#epoch');ylabel('Objective');
    title(['Convergence Speed and Runtime mat C #EXP ' num2str(eidx)]);
    set(gca,'fontsize',14);
    saveas(gca,[EXP.output_dir 'figure5.png']);
end






end