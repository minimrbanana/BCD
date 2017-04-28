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
KKT_A = zeros(3*EXP.n_loop,iters);
KKT_B = zeros(3*EXP.n_loop,iters);
KKT_C = zeros(3*EXP.n_loop,iters);
% to save objective, which is function value
OBJ_A = zeros(3*EXP.n_loop,iters);
OBJ_B = zeros(3*EXP.n_loop,iters);
OBJ_C = zeros(3*EXP.n_loop,iters);
%% loop to average the cenvergence
for loop=1:EXP.n_loop
    b = randn(EXP.d,1);
    Bb= b(perm);
    Cb = Bb(reod);
    %%
    % runtime of matrix A
    %
    profile on;
    [~, kkt1a] = CBCD_size1_mex_sparse(EXP.A, b, d, iters);
    [~, kkt2a] = CBCD_size2_ss(EXP.A, b, d, iters);
    [~, kkt3a] = CBCD_size3_ss(EXP.A, b, d, iters);
    p=profile('info');
    profile off;
    % save time for matrix A
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size1_mex_sparse')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(1,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(1,loop));
            fprintf('#epochs : %d \n',length(kkt1a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(2,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(2,loop));
            fprintf('#epochs : %d \n',length(kkt2a)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(3,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(3,loop));
            fprintf('#epochs : %d \n',length(kkt3a)-1);
        end
    end
    % save KKT condition, related to normal cone
    KKT_A(loop,1:length(kkt1a)) = kkt1a;
    KKT_A(loop+EXP.n_loop,1:length(kkt2a)) = kkt2a;
    KKT_A(loop+EXP.n_loop*2,1:length(kkt3a)) = kkt3a;
    % run again to get objective(function value)
    [~, obj1a] = CBCD_size1_fx(EXP.A, b, d, iters);
    [~, obj2a] = CBCD_size2_fx(EXP.A, b, d, iters);
    [~, obj3a] = CBCD_size3_fx(EXP.A, b, d, iters);
    % save objective, which is function value
    OBJ_A(loop,1:length(obj1a)) = obj1a;
    OBJ_A(loop+EXP.n_loop,1:length(obj2a)) = obj2a;
    OBJ_A(loop+EXP.n_loop*2,1:length(obj3a)) = obj3a;
    %%
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
            fprintf('#epochs : %d \n',length(y1b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(5,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(5,loop));
            fprintf('#epochs : %d \n',length(y2b)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(6,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(6,loop));
            fprintf('#epochs : %d \n',length(y3b)-1);
        end
    end
    % save objective
    KKT_B(loop,1:length(y1b)) = y1b;
    KKT_B(loop+EXP.n_loop,1:length(y2b)) = y2b;
    KKT_B(loop+EXP.n_loop*2,1:length(y3b)) = y3b;
    % run again to get objective(function value)
    [~, obj1b] = CBCD_size1_fx(B, Bb, d, iters);
    [~, obj2b] = CBCD_size2_fx(B, Bb, d, iters);
    [~, obj3b] = CBCD_size3_fx(B, Bb, d, iters);
    % save objective, which is function value
    OBJ_B(loop,1:length(obj1b)) = obj1b;
    OBJ_B(loop+EXP.n_loop,1:length(obj2b)) = obj2b;
    OBJ_B(loop+EXP.n_loop*2,1:length(obj3b)) = obj3b;
    %%
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
            fprintf('#epochs : %d \n',length(y1c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size2_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(8,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(8,loop));
            fprintf('#epochs : %d \n',length(y2c)-1);
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'CBCD_size3_ss')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(9,loop)= p.FunctionTable(i,1).TotalTime;
            fprintf('Runtime : %.4f seconds.   ',T(9,loop));
            fprintf('#epochs : %d \n',length(y3c)-1);
        end
    end
    % save objective
    KKT_C(loop,1:length(y1c)) = y1c;
    KKT_C(loop+EXP.n_loop,1:length(y2c)) = y2c;
    KKT_C(loop+EXP.n_loop*2,1:length(y3c)) = y3c;
    % run again to get objective(function value)
    [~, obj1c] = CBCD_size1_fx(C, Cb, d, iters);
    [~, obj2c] = CBCD_size2_fx(C, Cb, d, iters);
    [~, obj3c] = CBCD_size3_fx(C, Cb, d, iters);
    % save objective, which is function value
    OBJ_C(loop,1:length(obj1c)) = obj1c;
    OBJ_C(loop+EXP.n_loop,1:length(obj2c)) = obj2c;
    OBJ_C(loop+EXP.n_loop*2,1:length(obj3c)) = obj3c;
end
%% this part shows the mean of ratio of number of epochs 
% taken by size2/3 over size1, for matrix A,B and C
% and then remove the all-0 columns in KKT and OBJ
epoch1 = sum(KKT_A~=0,2);
epoch2 = sum(KKT_B~=0,2);
epoch3 = sum(KKT_C~=0,2);
% remove the experiments that reaches the max iters
index2remove = find(epoch1(1:EXP.n_loop)==iters);
% construct the index of size 123
vector2remove = [index2remove;index2remove+EXP.n_loop;index2remove+EXP.n_loop*2];
epoch1(vector2remove)=[];
epoch2(vector2remove)=[];
epoch3(vector2remove)=[];
KKT_A(vector2remove,:)=[];
KKT_B(vector2remove,:)=[];
KKT_C(vector2remove,:)=[];
OBJ_A(vector2remove,:)=[];
OBJ_B(vector2remove,:)=[];
OBJ_C(vector2remove,:)=[];
EXP.n_loop = size(epoch1,1)/3;
EXP.A_ep2_OVER_ep1 = mean(epoch1(1+EXP.n_loop:2*EXP.n_loop)./epoch1(1:EXP.n_loop));
EXP.A_ep3_OVER_ep1 = mean(epoch1(1+2*EXP.n_loop:3*EXP.n_loop)./epoch1(1:EXP.n_loop));
EXP.B_ep2_OVER_ep1 = mean(epoch2(1+EXP.n_loop:2*EXP.n_loop)./epoch2(1:EXP.n_loop));
EXP.B_ep3_OVER_ep1 = mean(epoch2(1+2*EXP.n_loop:3*EXP.n_loop)./epoch2(1:EXP.n_loop));
EXP.C_ep2_OVER_ep1 = mean(epoch3(1+EXP.n_loop:2*EXP.n_loop)./epoch3(1:EXP.n_loop));
EXP.C_ep3_OVER_ep1 = mean(epoch3(1+2*EXP.n_loop:3*EXP.n_loop)./epoch3(1:EXP.n_loop));
% to check whether reaches the max number of iterations or not
EXP.check_iter = max([epoch1;epoch2;epoch3]);

%% save the parameters to EXP
EXP.T = T;
EXP.KKT_A = KKT_A(:,1:max(epoch1));
EXP.KKT_B = KKT_B(:,1:max(epoch2));
EXP.KKT_C = KKT_C(:,1:max(epoch3));
EXP.OBJ_A = OBJ_A(:,1:max(epoch1));
EXP.OBJ_B = OBJ_B(:,1:max(epoch2));
EXP.OBJ_C = OBJ_C(:,1:max(epoch3));
EXP.B = B;
EXP.C = C;
% to save EXP
if EXP.save==1
    save([EXP.output_dir 'EXP.mat'],'EXP');
end
% if plot the convergence
if EXP.plot_convergence ==1
    plot4EXP(eidx,EXP);
end

end