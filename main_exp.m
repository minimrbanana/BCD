function EXP = main_exp(eidx)
% the main function for the experiments and save the results
%% clear and then load the experimental parameters
close all;
clear EXP;
EXP = exp_detail(eidx);
d = EXP.d;
iters = EXP.max_iter;
% do permutation on A, we get B, reordering on B we get C
perm = randperm(EXP.d);
B = EXP.A(perm,perm);
% reordering and estimate its run time
tstart = tic;
reod = symrcm(B);
C = B(reod,reod);
EXP.tRCM = toc(tstart);

% create variables for saving the results
% save runtime
T = zeros(9,EXP.n_loop);
% to save the KKT condition, related to normal cone
% 1,2,3 for matrix A, B, C
% first n_loop rows store size1
% second n_loop rows store size2
% third n_loop rows store size3
KKT_A = cell(3*EXP.n_loop,1);
KKT_B = cell(3*EXP.n_loop,1);
KKT_C = cell(3*EXP.n_loop,1);
% to save objective, which is function value
OBJ_A = cell(3*EXP.n_loop,1);
OBJ_B = cell(3*EXP.n_loop,1);
OBJ_C = cell(3*EXP.n_loop,1);
% to save the number of epochs
epoch1 = zeros(3*EXP.n_loop,1);
epoch2 = zeros(3*EXP.n_loop,1);
epoch3 = zeros(3*EXP.n_loop,1);
%% loop to average the cenvergence
for loop=1:EXP.n_loop
    b = randn(EXP.d,1);
    Bb= b(perm);
    Cb = Bb(reod);
    %%
    % runtime of matrix A
    %
    profile on;
    [~, kkt1a] = CBCD_size1_gc(EXP.A, b, d, iters,1E-13,0,1,0);
    [~, kkt2a] = CBCD_size2_gc(EXP.A, b, d, iters,1E-13,0,1,0);
    [~, kkt3a] = CBCD_size3_gc(EXP.A, b, d, iters,1E-13,0,1,0);
    p=profile('info');
    profile off;
    % save time for matrix A
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(1,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(1,loop));
            fprintf('#epochs : %d \n',length(kkt1a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(2,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(2,loop));
            fprintf('#epochs : %d \n',length(kkt2a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(3,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(3,loop));
            fprintf('#epochs : %d \n',length(kkt3a)-1);
        end
    end
    % save KKT condition, related to normal cone
    KKT_A{loop} = kkt1a;
    epoch1(loop) = length(kkt1a);
    KKT_A{loop+EXP.n_loop} = kkt2a;
    epoch1(loop+EXP.n_loop) = length(kkt2a);
    KKT_A{loop+EXP.n_loop*2} = kkt3a;
    epoch1(loop+EXP.n_loop*2) = length(kkt3a);
    % run again to get objective(function value)
    [~, obj1a] = CBCD_size1_fx(EXP.A, b, d, iters);
    [~, obj2a] = CBCD_size2_fx(EXP.A, b, d, iters);
    [~, obj3a] = CBCD_size3_fx(EXP.A, b, d, iters);
    % save objective, which is function value
    OBJ_A{loop} = obj1a;
    OBJ_A{loop+EXP.n_loop} = obj2a;
    OBJ_A{loop+EXP.n_loop*2} = obj3a;
    %%
    % runtime of matrix B
    %
    profile on;
    [~, kkt1b] = CBCD_size1_gc(B, Bb, d, iters,1E-13,0,1,0);
    [~, kkt2b] = CBCD_size2_gc(B, Bb, d, iters,1E-13,0,1,0);
    [~, kkt3b] = CBCD_size3_gc(B, Bb, d, iters,1E-13,0,1,0);
    p=profile('info');
    profile off;
    % save time for matrix B
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(4,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(4,loop));
            fprintf('#epochs : %d \n',length(kkt1b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(5,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(5,loop));
            fprintf('#epochs : %d \n',length(kkt2b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(6,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(6,loop));
            fprintf('#epochs : %d \n',length(kkt3b)-1);
        end
    end
    % save objective
    KKT_B{loop} = kkt1b;
    epoch2(loop) = length(kkt1b);
    KKT_B{loop+EXP.n_loop} = kkt2b;
    epoch2(loop+EXP.n_loop) = length(kkt2b);
    KKT_B{loop+EXP.n_loop*2} = kkt3b;
    epoch2(loop+EXP.n_loop*2) = length(kkt3b);
    % run again to get objective(function value)
    [~, obj1b] = CBCD_size1_fx(B, Bb, d, iters);
    [~, obj2b] = CBCD_size2_fx(B, Bb, d, iters);
    [~, obj3b] = CBCD_size3_fx(B, Bb, d, iters);
    % save objective, which is function value
    OBJ_B{loop} = obj1b;
    OBJ_B{loop+EXP.n_loop} = obj2b;
    OBJ_B{loop+EXP.n_loop*2} = obj3b;
    %%
    % runtime of matrix C
    %
    profile on;
    [~, kkt1c] = CBCD_size1_gc(C, Cb, d, iters,1E-13,0,1,0);
    [~, kkt2c] = CBCD_size2_gc(C, Cb, d, iters,1E-13,0,1,0);
    [~, kkt3c] = CBCD_size3_gc(C, Cb, d, iters,1E-13,0,1,0);
    p=profile('info');
    profile off;
    % save time for matrix C
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(7,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(7,loop));
            fprintf('#epochs : %d \n',length(kkt1c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(8,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(8,loop));
            fprintf('#epochs : %d \n',length(kkt2c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(9,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(9,loop));
            fprintf('#epochs : %d \n',length(kkt3c)-1);
        end
    end
    % save objective
    KKT_C{loop} = kkt1c;
    epoch3(loop) = length(kkt1c);
    KKT_C{loop+EXP.n_loop} = kkt2c;
    epoch3(loop+EXP.n_loop) = length(kkt2c);
    KKT_C{loop+EXP.n_loop*2} = kkt3c;
    epoch3(loop+EXP.n_loop*2) = length(kkt3a);
    % run again to get objective(function value)
    [~, obj1c] = CBCD_size1_fx(C, Cb, d, iters);
    [~, obj2c] = CBCD_size2_fx(C, Cb, d, iters);
    [~, obj3c] = CBCD_size3_fx(C, Cb, d, iters);
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

%% save the parameters to EXP
EXP.T = T;
EXP.KKT_A = KKT_A;
EXP.KKT_B = KKT_B;
EXP.KKT_C = KKT_C;
EXP.OBJ_A = OBJ_A;
EXP.OBJ_B = OBJ_B;
EXP.OBJ_C = OBJ_C;
EXP.B = B;
EXP.C = C;
EXP.epoch1 = epoch1;
EXP.epoch2 = epoch2;
EXP.epoch3 = epoch3;
% to save EXP
if EXP.save==1
    save([EXP.output_dir 'EXP.mat'],'EXP');
end
% if plot the convergence
if EXP.plot_convergence ==1
    plot4EXP(eidx,EXP);
end

end