%% percentage in the block and convergence
% plot the convergence of band matrix with 0% to 100% permutation
% as well as the |M|_2 of these matrices

% plot both 3-diagonal and 5,7,9-diagonal band matrix A 

% setting
rng(1);
addpath ../
%% parameters
% input
d = 300; 
lambda=1E-3;
iters=200000;
N=10000;
Repeat=10;
% output
epoch1=zeros(11,1);
epoch2=zeros(11,1);
epoch3=zeros(11,1);
cover1=zeros(11,1);
cover2=zeros(11,1);
cover3=zeros(11,1);
Mnorm1=zeros(11,1);
Mnorm2=zeros(11,1);
Mnorm3=zeros(11,1);
Meig1=zeros(11,1);
Meig2=zeros(11,1);
Meig3=zeros(11,1);
% M with 0%,30%,60% 100% permutation, for plotting ||M^{t}||_{2}
M1_0=0;
M1_30=0;
M1_60=0;
M1_100=0;
M2_0=0;
M2_30=0;
M2_60=0;
M2_100=0;
M3_0=0;
M3_30=0;
M3_60=0;
M3_100=0;
%% data 3-diagonal matrix
e1 = ones(d,1);
lambdaI=speye(d)*lambda;
A=spdiags([-e1,2*e1,-e1],[-1,0,1],d,d);
A(1,d)=-1;A(d,1)=-1;
A=A+lambdaI;
b=zeros(d,1);
nz=nnz(A);
% template
L = tril(ones(d,d),1);
L1=L;L2=L;L3=L;
for i=1:d-1
    L1(i,i+1)=0;
    if mod(i,2)==0
        L2(i,i+1)=0;
    end
    if mod(i,3)==0
        L3(i,i+1)=0;
    end
end

block1=spdiags(ones(d,1),0,d,d);

e0 = [1;0];
e1 = e0(:,ones(ceil(d/2),1));
e1 = reshape(e1 ,numel(e1),1);
block2 = spdiags([e1,[0;e1(1:end-1)]],[-1,1],d,d)+block1;

e0 = [1;1;0];
e1 = e0(:,ones(ceil(d/3),1));
e1 = reshape(e1 ,numel(e1),1);
e0 = [1;0;0];
e2 = e0(:,ones(ceil(d/3),1));
e2 = reshape(e2 ,numel(e2),1);
block3 = spdiags([e2,e1,[0;e1(1:end-1)],[0;0;e2(1:end-2)]],...
            [-2,-1,1,2],d,d)+block1;
% original matrix A
% solver
[cx01, cy1] = CBCD1(A, b, d, iters,1E-10,-100,100,1);
[cx02, cy2] = CBCD2(A, b, d, iters,1E-10,-100,100,1);
[cx03, cy3] = CBCD3(A, b, d, iters,1E-10,-100,100,1);
% save #epochs
epoch1(1)=size(cy1,1)*Repeat;
epoch2(1)=size(cy2,1)*Repeat;
epoch3(1)=size(cy3,1)*Repeat;
% save cover
cover1(1)=d/nz;
cover2(1)=nnz(A.*block2)/nz;
cover3(1)=nnz(A.*block3)/nz;
% M
I=eye(d);
AL1=A.*L1;
M1=I-AL1\A;Mnorm1(1)=norm(M1^N)*Repeat;
Meig1(1)=max(abs(eig(M1)))*Repeat;
AL2=A.*L2;
M2=I-AL2\A;Mnorm2(1)=norm(M2^N)*Repeat;
Meig2(1)=max(abs(eig(M2)))*Repeat;
AL3=A.*L3;
M3=I-AL3\A;Mnorm3(1)=norm(M3^N)*Repeat;
Meig3(1)=max(abs(eig(M3)))*Repeat;
M1_0=M1;
M2_0=M2;
M3_0=M3;

h3=figure(3);clf,suptitle('Structure of Matrix A');
set(gca,'fontsize',14);
%repeat 10 times
for j=1:Repeat
% 10 shuffle
for loop=1:10
    p0=randperm(round(d/10)*loop);
    p =[p0,round(d/10)*loop+1:d];
    A1=A(p,p);
    % solver
    [~, cy1] = CBCD1(A1, b, d, iters,1E-10,-100,100,1);
    [~, cy2] = CBCD2(A1, b, d, iters,1E-10,-100,100,1);
    [~, cy3] = CBCD3(A1, b, d, iters,1E-10,-100,100,1);
    % save #epochs
    epoch1(loop+1)=epoch1(loop+1)+size(cy1,1);
    epoch2(loop+1)=epoch2(loop+1)+size(cy2,1);
    epoch3(loop+1)=epoch3(loop+1)+size(cy3,1);
    % save cover
    cover1(loop+1)=d/nz;
    cover2(loop+1)=nnz(A1.*block2)/nz;
    cover3(loop+1)=nnz(A1.*block3)/nz;
    % M
    I=eye(d);
    AL1=A1.*L1;
    M1=I-AL1\A1;
    Mnorm1(loop+1)=Mnorm1(loop+1)+norm(M1^N);
    Meig1(loop+1) =Meig1(loop+1)+max(abs(eig(M1)));
    AL2=A1.*L2;
    M2=I-AL2\A1;
    Mnorm2(loop+1)=Mnorm2(loop+1)+norm(M2^N);
    Meig2(loop+1) =Meig2(loop+1)+max(abs(eig(M2)));
    AL3=A1.*L3;
    M3=I-AL3\A1;
    Mnorm3(loop+1)=Mnorm3(loop+1)+norm(M3^N);
    Meig3(loop+1) =Meig3(loop+1)+max(abs(eig(M3)));
    if loop<10
        subplot(3,3,loop),
        spy(A1);
    end
    if M1_30==0
        if loop==3
            M1_30=M1;
            M2_30=M2;
            M3_30=M3;
        end
    end
    if M1_60==0
        if loop==6
            M1_60=M1;
            M2_60=M2;
            M3_60=M3;
        end
    end
    if M1_100==0
        if loop==10
            M1_100=M1;
            M2_100=M2;
            M3_100=M3;
        end
    end
end

end

save('/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band3.mat',...
    'cover1','cover2','cover3','epoch1','epoch2','epoch3','Meig1','Meig2','Meig3','Repeat');

h2=figure(2);clf,
set(gca,'fontsize',14);
subplot(2,4,1),
cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
plot(cover1,epoch1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD');
axis([0 1 0.5E4 2.6E4]);
set(gca,'FontSize',12);

% subplot(3,4,5),
% % plot(cover1,Mnorm1,'ro-');hold on;
% % plot(cover2,Mnorm2,'go-');hold on;
% % plot(cover3,Mnorm3,'bo-');hold on;
% semilogy(cover1,Mnorm1/Repeat,'ro-','LineWidth',1.5);hold on;
% semilogy(cover2,Mnorm2/Repeat,'go-','LineWidth',1.5);hold on;
% semilogy(cover3,Mnorm3/Repeat,'bo-','LineWidth',1.5);hold on;
% legend('CBCD1','CBCD2','CBCD3','Location','southeast');
% grid on;
% xlabel('Permutation rate');
% ylabel('||M^{10000}||_{2}');
% title('Matrix Norm of M^{10000}');

subplot(2,4,5),
plot(cover1,Meig1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
set(gca,'FontSize',12);

%% 5 band matrix
N=20000;
% output
epoch1=zeros(11,1);
epoch2=zeros(11,1);
epoch3=zeros(11,1);
cover1=zeros(11,1);
cover2=zeros(11,1);
cover3=zeros(11,1);
Mnorm1=zeros(11,1);
Mnorm2=zeros(11,1);
Mnorm3=zeros(11,1);
Meig1=zeros(11,1);
Meig2=zeros(11,1);
Meig3=zeros(11,1);
e1 = ones(d,1);
lambdaI=speye(d)*lambda;
A=spdiags([-e1,-e1,4*e1,-e1,-e1],[-2,-1,0,1,2],d,d);
A(1,d)=-1;A(d,1)=-1;
A=A+lambdaI;
b=zeros(d,1);
nz=nnz(A);
% template
L = tril(ones(d,d),1);
L1=L;L2=L;L3=L;
for i=1:d-1
    L1(i,i+1)=0;
    if mod(i,2)==0
        L2(i,i+1)=0;
    end
    if mod(i,3)==0
        L3(i,i+1)=0;
    end
end

block1=spdiags(ones(d,1),0,d,d);

e0 = [1;0];
e1 = e0(:,ones(ceil(d/2),1));
e1 = reshape(e1 ,numel(e1),1);
block2 = spdiags([e1,[0;e1(1:end-1)]],[-1,1],d,d)+block1;

e0 = [1;1;0];
e1 = e0(:,ones(ceil(d/3),1));
e1 = reshape(e1 ,numel(e1),1);
e0 = [1;0;0];
e2 = e0(:,ones(ceil(d/3),1));
e2 = reshape(e2 ,numel(e2),1);
block3 = spdiags([e2,e1,[0;e1(1:end-1)],[0;0;e2(1:end-2)]],...
            [-2,-1,1,2],d,d)+block1;
% original matrix A
% solver
[~, cy1] = CBCD1(A, b, d, iters,1E-10,-100,100,1);
[~, cy2] = CBCD2(A, b, d, iters,1E-10,-100,100,1);
[~, cy3] = CBCD3(A, b, d, iters,1E-10,-100,100,1);
% save #epochs
epoch1(1)=size(cy1,1)*Repeat;
epoch2(1)=size(cy2,1)*Repeat;
epoch3(1)=size(cy3,1)*Repeat;
% save cover
cover1(1)=d/nz;
cover2(1)=nnz(A.*block2)/nz;
cover3(1)=nnz(A.*block3)/nz;
% M
I=eye(d);
AL1=A.*L1;
M1=I-AL1\A;Mnorm1(1)=norm(M1^N)*Repeat;
Meig1(1)=max(abs(eig(M1)))*Repeat;
AL2=A.*L2;
M2=I-AL2\A;Mnorm2(1)=norm(M2^N)*Repeat;
Meig2(1)=max(abs(eig(M2)))*Repeat;
AL3=A.*L3;
M3=I-AL3\A;Mnorm3(1)=norm(M3^N)*Repeat;
Meig3(1)=max(abs(eig(M3)))*Repeat;

% repeat 10 times
for j=1:Repeat
% 10 shuffle
for loop=1:10
    p0=randperm(round(d/10)*loop);
    p =[p0,round(d/10)*loop+1:d];
    A1=A(p,p);
    % solver
    [cx01, cy1] = CBCD1(A1, b, d, iters,1E-10,-100,100,1);
    [cx02, cy2] = CBCD2(A1, b, d, iters,1E-10,-100,100,1);
    [cx03, cy3] = CBCD3(A1, b, d, iters,1E-10,-100,100,1);
    % save #epochs
    epoch1(loop+1)=epoch1(loop+1)+size(cy1,1);
    epoch2(loop+1)=epoch2(loop+1)+size(cy2,1);
    epoch3(loop+1)=epoch3(loop+1)+size(cy3,1);
    % save cover
    cover1(loop+1)=d/nz;
    cover2(loop+1)=nnz(A1.*block2)/nz;
    cover3(loop+1)=nnz(A1.*block3)/nz;
    % M
    I=eye(d);
    AL1=A1.*L1;
    M1=I-AL1\A1;
    Mnorm1(loop+1)=Mnorm1(loop+1)+norm(M1^N);
    Meig1(loop+1) =Meig1(loop+1)+max(abs(eig(M1)));
    AL2=A1.*L2;
    M2=I-AL2\A1;
    Mnorm2(loop+1)=Mnorm2(loop+1)+norm(M2^N);
    Meig2(loop+1) =Meig2(loop+1)+max(abs(eig(M2)));
    AL3=A1.*L3;
    M3=I-AL3\A1;
    Mnorm3(loop+1)=Mnorm3(loop+1)+norm(M3^N);
    Meig3(loop+1) =Meig3(loop+1)+max(abs(eig(M3)));
end

end

save('/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band5.mat',...
    'cover1','cover2','cover3','epoch1','epoch2','epoch3','Meig1','Meig2','Meig3','Repeat');

cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
subplot(2,4,2),
plot(cover1,epoch1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD');
axis([0 1 0.5E4 2.6E4]);
set(gca,'FontSize',12);
% subplot(3,4,6),
% % plot(cover1,Mnorm1,'ro-');hold on;
% % plot(cover2,Mnorm2,'go-');hold on;
% % plot(cover3,Mnorm3,'bo-');hold on;
% semilogy(cover1,Mnorm1/Repeat,'ro-','LineWidth',1.5);hold on;
% semilogy(cover2,Mnorm2/Repeat,'go-','LineWidth',1.5);hold on;
% semilogy(cover3,Mnorm3/Repeat,'bo-','LineWidth',1.5);hold on;
% legend('CBCD1','CBCD2','CBCD3','Location','southeast');
% grid on;
% xlabel('Permutation rate');
% ylabel('||M^{20000}||_{2}');
% title('Matrix Norm of M^{20000}');

subplot(2,4,6),
plot(cover1,Meig1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
set(gca,'FontSize',12);
%% 7 band matrix
N=20000;
% output
epoch1=zeros(11,1);
epoch2=zeros(11,1);
epoch3=zeros(11,1);
cover1=zeros(11,1);
cover2=zeros(11,1);
cover3=zeros(11,1);
Mnorm1=zeros(11,1);
Mnorm2=zeros(11,1);
Mnorm3=zeros(11,1);
Meig1=zeros(11,1);
Meig2=zeros(11,1);
Meig3=zeros(11,1);
e1 = ones(d,1);
lambdaI=speye(d)*lambda;
A=spdiags([-e1,-e1,-e1,6*e1,-e1,-e1,-e1],[-3,-2,-1,0,1,2,3],d,d);
A(1,d)=-1;A(d,1)=-1;
A=A+lambdaI;
b=zeros(d,1);
nz=nnz(A);
% template
L = tril(ones(d,d),1);
L1=L;L2=L;L3=L;
for i=1:d-1
    L1(i,i+1)=0;
    if mod(i,2)==0
        L2(i,i+1)=0;
    end
    if mod(i,3)==0
        L3(i,i+1)=0;
    end
end

block1=spdiags(ones(d,1),0,d,d);

e0 = [1;0];
e1 = e0(:,ones(ceil(d/2),1));
e1 = reshape(e1 ,numel(e1),1);
block2 = spdiags([e1,[0;e1(1:end-1)]],[-1,1],d,d)+block1;

e0 = [1;1;0];
e1 = e0(:,ones(ceil(d/3),1));
e1 = reshape(e1 ,numel(e1),1);
e0 = [1;0;0];
e2 = e0(:,ones(ceil(d/3),1));
e2 = reshape(e2 ,numel(e2),1);
block3 = spdiags([e2,e1,[0;e1(1:end-1)],[0;0;e2(1:end-2)]],...
            [-2,-1,1,2],d,d)+block1;
% original matrix A
% solver
[~, cy1] = CBCD1(A, b, d, iters,1E-10,-100,100,1);
[~, cy2] = CBCD2(A, b, d, iters,1E-10,-100,100,1);
[~, cy3] = CBCD3(A, b, d, iters,1E-10,-100,100,1);
% save #epochs
epoch1(1)=size(cy1,1)*Repeat;
epoch2(1)=size(cy2,1)*Repeat;
epoch3(1)=size(cy3,1)*Repeat;
% save cover
cover1(1)=d/nz;
cover2(1)=nnz(A.*block2)/nz;
cover3(1)=nnz(A.*block3)/nz;
% M
I=eye(d);
AL1=A.*L1;
M1=I-AL1\A;Mnorm1(1)=norm(M1^N)*Repeat;
Meig1(1)=max(abs(eig(M1)))*Repeat;
AL2=A.*L2;
M2=I-AL2\A;Mnorm2(1)=norm(M2^N)*Repeat;
Meig2(1)=max(abs(eig(M2)))*Repeat;
AL3=A.*L3;
M3=I-AL3\A;Mnorm3(1)=norm(M3^N)*Repeat;
Meig3(1)=max(abs(eig(M3)))*Repeat;

% repeat 10 times
for j=1:Repeat
% 10 shuffle
for loop=1:10
    p0=randperm(round(d/10)*loop);
    p =[p0,round(d/10)*loop+1:d];
    A1=A(p,p);
    % solver
    [cx01, cy1] = CBCD1(A1, b, d, iters,1E-10,-100,100,1);
    [cx02, cy2] = CBCD2(A1, b, d, iters,1E-10,-100,100,1);
    [cx03, cy3] = CBCD3(A1, b, d, iters,1E-10,-100,100,1);
    % save #epochs
    epoch1(loop+1)=epoch1(loop+1)+size(cy1,1);
    epoch2(loop+1)=epoch2(loop+1)+size(cy2,1);
    epoch3(loop+1)=epoch3(loop+1)+size(cy3,1);
    % save cover
    cover1(loop+1)=d/nz;
    cover2(loop+1)=nnz(A1.*block2)/nz;
    cover3(loop+1)=nnz(A1.*block3)/nz;
    % M
    I=eye(d);
    AL1=A1.*L1;
    M1=I-AL1\A1;
    Mnorm1(loop+1)=Mnorm1(loop+1)+norm(M1^N);
    Meig1(loop+1) =Meig1(loop+1)+max(abs(eig(M1)));
    AL2=A1.*L2;
    M2=I-AL2\A1;
    Mnorm2(loop+1)=Mnorm2(loop+1)+norm(M2^N);
    Meig2(loop+1) =Meig2(loop+1)+max(abs(eig(M2)));
    AL3=A1.*L3;
    M3=I-AL3\A1;
    Mnorm3(loop+1)=Mnorm3(loop+1)+norm(M3^N);
    Meig3(loop+1) =Meig3(loop+1)+max(abs(eig(M3)));
end

end

save('/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band7.mat',...
    'cover1','cover2','cover3','epoch1','epoch2','epoch3','Meig1','Meig2','Meig3','Repeat');

cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
subplot(2,4,3),
plot(cover1,epoch1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD');
axis([0 1 0.5E4 2.6E4]);
set(gca,'FontSize',12);
% subplot(3,4,7),
% % plot(cover1,Mnorm1,'ro-');hold on;
% % plot(cover2,Mnorm2,'go-');hold on;
% % plot(cover3,Mnorm3,'bo-');hold on;
% semilogy(cover1,Mnorm1/Repeat,'ro-','LineWidth',1.5);hold on;
% semilogy(cover2,Mnorm2/Repeat,'go-','LineWidth',1.5);hold on;
% semilogy(cover3,Mnorm3/Repeat,'bo-','LineWidth',1.5);hold on;
% legend('CBCD1','CBCD2','CBCD3','Location','southeast');
% grid on;
% xlabel('Permutation rate');
% ylabel('||M^{20000}||_{2}');
% title('Matrix Norm of M^{20000}');

subplot(2,4,7),
plot(cover1,Meig1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
set(gca,'FontSize',12);

%% 9 band matrix
N=20000;
% output
epoch1=zeros(11,1);
epoch2=zeros(11,1);
epoch3=zeros(11,1);
cover1=zeros(11,1);
cover2=zeros(11,1);
cover3=zeros(11,1);
Mnorm1=zeros(11,1);
Mnorm2=zeros(11,1);
Mnorm3=zeros(11,1);
Meig1=zeros(11,1);
Meig2=zeros(11,1);
Meig3=zeros(11,1);
e1 = ones(d,1);
lambdaI=speye(d)*lambda;
%A=spdiags([-e1,-e1,4*e1,-e1,-e1],[-2,-1,0,1,2],d,d);
A=spdiags([-e1,-e1,-e1,-e1,8*e1,-e1,-e1,-e1,-e1],[-4,-3,-2,-1,0,1,2,3,4],d,d);
A(1,d)=-1;A(d,1)=-1;
A=A+lambdaI;
b=zeros(d,1);
nz=nnz(A);
% template
L = tril(ones(d,d),1);
L1=L;L2=L;L3=L;
for i=1:d-1
    L1(i,i+1)=0;
    if mod(i,2)==0
        L2(i,i+1)=0;
    end
    if mod(i,3)==0
        L3(i,i+1)=0;
    end
end

block1=spdiags(ones(d,1),0,d,d);

e0 = [1;0];
e1 = e0(:,ones(ceil(d/2),1));
e1 = reshape(e1 ,numel(e1),1);
block2 = spdiags([e1,[0;e1(1:end-1)]],[-1,1],d,d)+block1;

e0 = [1;1;0];
e1 = e0(:,ones(ceil(d/3),1));
e1 = reshape(e1 ,numel(e1),1);
e0 = [1;0;0];
e2 = e0(:,ones(ceil(d/3),1));
e2 = reshape(e2 ,numel(e2),1);
block3 = spdiags([e2,e1,[0;e1(1:end-1)],[0;0;e2(1:end-2)]],...
            [-2,-1,1,2],d,d)+block1;
% original matrix A
% solver
[~, cy1] = CBCD1(A, b, d, iters,1E-10,-100,100,1);
[~, cy2] = CBCD2(A, b, d, iters,1E-10,-100,100,1);
[~, cy3] = CBCD3(A, b, d, iters,1E-10,-100,100,1);
% save #epochs
epoch1(1)=size(cy1,1)*Repeat;
epoch2(1)=size(cy2,1)*Repeat;
epoch3(1)=size(cy3,1)*Repeat;
% save cover
cover1(1)=d/nz;
cover2(1)=nnz(A.*block2)/nz;
cover3(1)=nnz(A.*block3)/nz;
% M
I=eye(d);
AL1=A.*L1;
M1=I-AL1\A;Mnorm1(1)=norm(M1^N)*Repeat;
Meig1(1)=max(abs(eig(M1)))*Repeat;
AL2=A.*L2;
M2=I-AL2\A;Mnorm2(1)=norm(M2^N)*Repeat;
Meig2(1)=max(abs(eig(M2)))*Repeat;
AL3=A.*L3;
M3=I-AL3\A;Mnorm3(1)=norm(M3^N)*Repeat;
Meig3(1)=max(abs(eig(M3)))*Repeat;

% repeat 10 times
for j=1:Repeat
% 10 shuffle
for loop=1:10
    p0=randperm(round(d/10)*loop);
    p =[p0,round(d/10)*loop+1:d];
    A1=A(p,p);
    % solver
    [cx01, cy1] = CBCD1(A1, b, d, iters,1E-10,-100,100,1);
    [cx02, cy2] = CBCD2(A1, b, d, iters,1E-10,-100,100,1);
    [cx03, cy3] = CBCD3(A1, b, d, iters,1E-10,-100,100,1);
    % save #epochs
    epoch1(loop+1)=epoch1(loop+1)+size(cy1,1);
    epoch2(loop+1)=epoch2(loop+1)+size(cy2,1);
    epoch3(loop+1)=epoch3(loop+1)+size(cy3,1);
    % save cover
    cover1(loop+1)=d/nz;
    cover2(loop+1)=nnz(A1.*block2)/nz;
    cover3(loop+1)=nnz(A1.*block3)/nz;
    % M
    I=eye(d);
    AL1=A1.*L1;
    M1=I-AL1\A1;
    Mnorm1(loop+1)=Mnorm1(loop+1)+norm(M1^N);
    Meig1(loop+1) =Meig1(loop+1)+max(abs(eig(M1)));
    AL2=A1.*L2;
    M2=I-AL2\A1;
    Mnorm2(loop+1)=Mnorm2(loop+1)+norm(M2^N);
    Meig2(loop+1) =Meig2(loop+1)+max(abs(eig(M2)));
    AL3=A1.*L3;
    M3=I-AL3\A1;
    Mnorm3(loop+1)=Mnorm3(loop+1)+norm(M3^N);
    Meig3(loop+1) =Meig3(loop+1)+max(abs(eig(M3)));
end

end

save('/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band9.mat',...
    'cover1','cover2','cover3','epoch1','epoch2','epoch3','Meig1','Meig2','Meig3','Repeat');


cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
subplot(2,4,4),
plot(cover1,epoch1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD');
axis([0 1 0.5E4 2.6E4]);
set(gca,'FontSize',12);
% subplot(3,4,8),
% % plot(cover1,Mnorm1,'ro-');hold on;
% % plot(cover2,Mnorm2,'go-');hold on;
% % plot(cover3,Mnorm3,'bo-');hold on;
% semilogy(cover1,Mnorm1/Repeat,'ro-','LineWidth',1.5);hold on;
% semilogy(cover2,Mnorm2/Repeat,'go-','LineWidth',1.5);hold on;
% semilogy(cover3,Mnorm3/Repeat,'bo-','LineWidth',1.5);hold on;
% legend('CBCD1','CBCD2','CBCD3','Location','southeast');
% grid on;
% xlabel('Permutation rate');
% ylabel('||M^{20000}||_{2}');
% title('Matrix Norm of M^{20000}');


subplot(2,4,8),
plot(cover1,Meig1/Repeat,'ro-','LineWidth',1.5);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',1.5);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',1.5);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
set(gca,'FontSize',12);

