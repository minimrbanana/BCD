function EXP = main_exp(eidx)
% the main function for the experiments and save the results
%% clear and then load the experimental parameters
close all;
clear EXP;
EXP = exp_detail(eidx);
% load the parameters from EXP
d = EXP.d;
iters = EXP.max_iter;
l=EXP.lower;
u=EXP.upper;
x0=EXP.init;
alpha =EXP.alpha;
pre = EXP.precision;
% if EXP.A exists, which means A is not random
% do permutation on A, we get B, reordering on B we get C
if isfield(EXP,'A')
    perm = randperm(EXP.d);
    B = EXP.A(perm,perm);
    % reordering and estimate its run time
    tstart = tic;
    reod = symrcm(B);
    C = B(reod,reod);
    tRCM = toc(tstart);
    runMyMatA=0;
else
    tRCM=zeros(EXP.n_loop,1);
    runMyMatA=1;
end
% create variables for saving the results
% save runtime
T_c = zeros(9,EXP.n_loop);
T_r = zeros(9,EXP.n_loop);
% to save the KKT condition, related to normal cone
% 1,2,3 for matrix A, B, C
% first n_loop rows store size1
% second n_loop rows store size2
% third n_loop rows store size3
KKT_A = cell(3*EXP.n_loop,1);
KKT_B = cell(3*EXP.n_loop,1);
KKT_C = cell(3*EXP.n_loop,1);
% to save KKT condition of RBCD
KKTr_A = cell(3*EXP.n_loop,1);
KKTr_B = cell(3*EXP.n_loop,1);
KKTr_C = cell(3*EXP.n_loop,1);
% to save objective, which is function value
OBJ_A = cell(3*EXP.n_loop,1);
OBJ_B = cell(3*EXP.n_loop,1);
OBJ_C = cell(3*EXP.n_loop,1);
% to save the number of epochs
epoch1 = zeros(3*EXP.n_loop,1);
epoch2 = zeros(3*EXP.n_loop,1);
epoch3 = zeros(3*EXP.n_loop,1);
% to save the number of epochs of RBCD
epoch1r = zeros(3*EXP.n_loop,1);
epoch2r = zeros(3*EXP.n_loop,1);
epoch3r = zeros(3*EXP.n_loop,1);
%% loop to average the cenvergence
for loop=1:EXP.n_loop
    if runMyMatA==1
        EXP.A = MyMatA(eidx,EXP.d);
        perm = randperm(EXP.d);
        B = EXP.A(perm,perm);
        % reordering and estimate its run time
        tstart = tic;
        reod = symrcm(B);
        C = B(reod,reod);
        tRCM(loop) = toc(tstart);
    end
    b = randn(EXP.d,1);
    Bb= b(perm);
    Cb = Bb(reod);
    %%
    % runtime of matrix A
    %
    profile on;
    [~, kkt1a]  = CBCD_size1_gc(EXP.A, b, d, iters,pre,l,u,x0);
    [~, kkt2a]  = CBCD_size2_gc(EXP.A, b, d, iters,pre,l,u,x0);
    [~, kkt3a]  = CBCD_size3_gc(EXP.A, b, d, iters,pre,l,u,x0);
    [~, kktr1a] = RBCD_size1_gc_u(EXP.A, b, d, iters,pre,l,u,x0,alpha);
    [~, kktr2a] = RBCD2(EXP.A, b, d, iters,pre,l,u,x0,alpha);
    [~, kktr3a] = RBCD3(EXP.A, b, d, iters,pre,l,u,x0,alpha);
    if length(kkt1a)>iters-2
        save('A.mat','EXP','b');
    end
    p=profile('info');
    profile off;
    % save time for matrix A
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(1,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(1,loop));
            fprintf('#epochs : %d \n',length(kkt1a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(2,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(2,loop));
            fprintf('#epochs : %d \n',length(kkt2a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(3,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(3,loop));
            fprintf('#epochs : %d \n',length(kkt3a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size1_gc_u')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(1,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(1,loop));
            fprintf('#epochs : %d \n',length(kktr1a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD2')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(2,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(2,loop));
            fprintf('#epochs : %d \n',length(kktr2a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD3')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(3,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(3,loop));
            fprintf('#epochs : %d \n',length(kktr3a)-1);
        end
    end
    % save KKT condition, related to normal cone
    KKT_A{loop} = kkt1a;
    epoch1(loop) = length(kkt1a);
    KKT_A{loop+EXP.n_loop} = kkt2a;
    epoch1(loop+EXP.n_loop) = length(kkt2a);
    KKT_A{loop+EXP.n_loop*2} = kkt3a;
    epoch1(loop+EXP.n_loop*2) = length(kkt3a);
    % save KKT condition of RBCD
    KKTr_A{loop} = kktr1a;
    epoch1r(loop) = length(kktr1a);
    KKTr_A{loop+EXP.n_loop} = kktr2a;
    epoch1r(loop+EXP.n_loop) = length(kktr2a);
    KKTr_A{loop+EXP.n_loop*2} = kktr3a;
    epoch1r(loop+EXP.n_loop*2) = length(kktr3a);
    % run again to get objective(function value)
    [~, obj1a] = CBCD_size1_fx(EXP.A, b, d, iters,pre,l,u,x0);
    [~, obj2a] = CBCD_size2_fx(EXP.A, b, d, iters,pre,l,u,x0);
    [~, obj3a] = CBCD_size3_fx(EXP.A, b, d, iters,pre,l,u,x0);
    % save objective, which is function value
    OBJ_A{loop} = obj1a;
    OBJ_A{loop+EXP.n_loop} = obj2a;
    OBJ_A{loop+EXP.n_loop*2} = obj3a;
    %%
    % runtime of matrix B
    %
    profile on;
    [~, kkt1b]  = CBCD_size1_gc(B, Bb, d, iters,pre,l,u,x0);
    [~, kkt2b]  = CBCD_size2_gc(B, Bb, d, iters,pre,l,u,x0);
    [~, kkt3b]  = CBCD_size3_gc(B, Bb, d, iters,pre,l,u,x0);
    [~, kktr1b] = RBCD_size1_gc_u(B, Bb, d, iters,pre,l,u,x0,alpha);
    [~, kktr2b] = RBCD2(B, Bb, d, iters,pre,l,u,x0,alpha);
    [~, kktr3b] = RBCD3(B, Bb, d, iters,pre,l,u,x0,alpha);
    p=profile('info');
    profile off;
    % save time for matrix B
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(4,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(4,loop));
            fprintf('#epochs : %d \n',length(kkt1b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(5,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(5,loop));
            fprintf('#epochs : %d \n',length(kkt2b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(6,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(6,loop));
            fprintf('#epochs : %d \n',length(kkt3b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size1_gc_u')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(4,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(4,loop));
            fprintf('#epochs : %d \n',length(kktr1b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD2')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(5,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(5,loop));
            fprintf('#epochs : %d \n',length(kktr2b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD3')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(6,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(6,loop));
            fprintf('#epochs : %d \n',length(kktr3b)-1);
        end
    end
    % save objective
    KKT_B{loop} = kkt1b;
    epoch2(loop) = length(kkt1b);
    KKT_B{loop+EXP.n_loop} = kkt2b;
    epoch2(loop+EXP.n_loop) = length(kkt2b);
    KKT_B{loop+EXP.n_loop*2} = kkt3b;
    epoch2(loop+EXP.n_loop*2) = length(kkt3b);
    % save objective of RBCD
    KKTr_B{loop} = kktr1b;
    epoch2r(loop) = length(kktr1b);
    KKTr_B{loop+EXP.n_loop} = kktr2b;
    epoch2r(loop+EXP.n_loop) = length(kktr2b);
    KKTr_B{loop+EXP.n_loop*2} = kktr3b;
    epoch2r(loop+EXP.n_loop*2) = length(kktr3b);
    % run again to get objective(function value)
    [~, obj1b] = CBCD_size1_fx(B, Bb, d, iters,pre,l,u,x0);
    [~, obj2b] = CBCD_size2_fx(B, Bb, d, iters,pre,l,u,x0);
    [~, obj3b] = CBCD_size3_fx(B, Bb, d, iters,pre,l,u,x0);
    % save objective, which is function value
    OBJ_B{loop} = obj1b;
    OBJ_B{loop+EXP.n_loop} = obj2b;
    OBJ_B{loop+EXP.n_loop*2} = obj3b;
    %%
    % runtime of matrix C
    %
    profile on;
    [~, kkt1c]  = CBCD_size1_gc(C, Cb, d, iters,pre,l,u,x0);
    [~, kkt2c]  = CBCD_size2_gc(C, Cb, d, iters,pre,l,u,x0);
    [~, kkt3c]  = CBCD_size3_gc(C, Cb, d, iters,pre,l,u,x0);
    [~, kktr1c] = RBCD_size1_gc_u(C, Cb, d, iters,pre,l,u,x0,alpha);
    [~, kktr2c] = RBCD2(C, Cb, d, iters,pre,l,u,x0,alpha);
    [~, kktr3c] = RBCD3(C, Cb, d, iters,pre,l,u,x0,alpha);
    p=profile('info');
    profile off;
    % save time for matrix C
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(7,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(7,loop));
            fprintf('#epochs : %d \n',length(kkt1c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(8,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(8,loop));
            fprintf('#epochs : %d \n',length(kkt2c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_c(9,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_c(9,loop));
            fprintf('#epochs : %d \n',length(kkt3c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size1_gc_u')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(7,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(7,loop));
            fprintf('#epochs : %d \n',length(kktr1c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD2')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(8,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(8,loop));
            fprintf('#epochs : %d \n',length(kktr2c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD3')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T_r(9,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T_r(9,loop));
            fprintf('#epochs : %d \n',length(kktr3c)-1);
        end
    end
    % save objective
    KKT_C{loop} = kkt1c;
    epoch3(loop) = length(kkt1c);
    KKT_C{loop+EXP.n_loop} = kkt2c;
    epoch3(loop+EXP.n_loop) = length(kkt2c);
    KKT_C{loop+EXP.n_loop*2} = kkt3c;
    epoch3(loop+EXP.n_loop*2) = length(kkt3c);
    % save objective of RBCD
    KKTr_C{loop} = kktr1c;
    epoch3r(loop) = length(kktr1c);
    KKTr_C{loop+EXP.n_loop} = kktr2c;
    epoch3r(loop+EXP.n_loop) = length(kktr2c);
    KKTr_C{loop+EXP.n_loop*2} = kktr3c;
    epoch3r(loop+EXP.n_loop*2) = length(kktr3c);
    % run again to get objective(function value)
    [~, obj1c] = CBCD_size1_fx(C, Cb, d, iters,pre,l,u,x0);
    [~, obj2c] = CBCD_size2_fx(C, Cb, d, iters,pre,l,u,x0);
    [~, obj3c] = CBCD_size3_fx(C, Cb, d, iters,pre,l,u,x0);
    % save objective, which is function value
    OBJ_C{loop} = obj1c;
    OBJ_C{loop+EXP.n_loop} = obj2c;
    OBJ_C{loop+EXP.n_loop*2} = obj3c;
end
%% this part shows the mean of ratio of number of epochs 
% taken by size2/3 over size1, for matrix A,B and C
% we don't take the experiments that reaches the max iters
% when calculating the #2/#1 and #3/#1.
% However they are counted when plotting in function plot4EXP()
index2remove = find(epoch1(1:EXP.n_loop)==iters);
% construct the index of size 123
vector2remove = [index2remove;index2remove+EXP.n_loop;index2remove+EXP.n_loop*2];
temp_epoch1 = epoch1;
temp_epoch2 = epoch2;
temp_epoch3 = epoch3;
temp_epoch1(vector2remove)=[];
temp_epoch2(vector2remove)=[];
temp_epoch3(vector2remove)=[];
temp_exp_loop=EXP.n_loop;
EXP.n_loop = size(temp_epoch1,1)/3;
EXP.A_ep2_OVER_ep1 = mean(temp_epoch1(1+EXP.n_loop:2*EXP.n_loop)./...
    temp_epoch1(1:EXP.n_loop));
EXP.A_ep3_OVER_ep1 = mean(temp_epoch1(1+2*EXP.n_loop:3*EXP.n_loop)./...
    temp_epoch1(1:EXP.n_loop));
EXP.B_ep2_OVER_ep1 = mean(temp_epoch2(1+EXP.n_loop:2*EXP.n_loop)./...
    temp_epoch2(1:EXP.n_loop));
EXP.B_ep3_OVER_ep1 = mean(temp_epoch2(1+2*EXP.n_loop:3*EXP.n_loop)./...
    temp_epoch2(1:EXP.n_loop));
EXP.C_ep2_OVER_ep1 = mean(temp_epoch3(1+EXP.n_loop:2*EXP.n_loop)./...
    temp_epoch3(1:EXP.n_loop));
EXP.C_ep3_OVER_ep1 = mean(temp_epoch3(1+2*EXP.n_loop:3*EXP.n_loop)./...
    temp_epoch3(1:EXP.n_loop));
% to check whether reaches the max number of iterations or not
EXP.check_iter = max([temp_epoch1;temp_epoch2;temp_epoch3]);
% set the number of loops back to original
EXP.n_loop=temp_exp_loop;
%% save the parameters to EXP
EXP.T_c = T_c;
EXP.T_r = T_r;
EXP.KKT_A = KKT_A;
EXP.KKT_B = KKT_B;
EXP.KKT_C = KKT_C;
EXP.KKTr_A = KKTr_A;
EXP.KKTr_B = KKTr_B;
EXP.KKTr_C = KKTr_C;
EXP.OBJ_A = OBJ_A;
EXP.OBJ_B = OBJ_B;
EXP.OBJ_C = OBJ_C;
EXP.B = B;
EXP.C = C;
EXP.epoch1 = epoch1;
EXP.epoch2 = epoch2;
EXP.epoch3 = epoch3;
EXP.epoch1r = epoch1r;
EXP.epoch2r = epoch2r;
EXP.epoch3r = epoch3r;
EXP.tRCM=mean(tRCM);
% to save EXP
if EXP.save==1
    save([EXP.output_dir 'EXP.mat'],'EXP');
end
% if plot the convergence
if EXP.plot_convergence ==1
    plot4EXP(eidx,EXP);
end

end
