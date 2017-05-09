%% BDC with block size 1 & 2 & 3

%% input
%rng(1);
% dimension and constraints
d = 50;
%lower = zeros(d,1);
%upper = ones(d,1);
% A and b

e1 = ones(d,1);
A = spdiags([-e1,-e1],[-1,1],d,d);
diagonal = -sum(A);
diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
A = spdiags(diagonal',0,A);

% c = 0.8;
% A = ones(d,d)*c;
% A = A + diag(ones(d,1)*(1-c));
% A = sparse(A);
%A = sparse(A);
% b = zeros(d,1);
%rng(1);
b = randn(d,1);
%% solution by dense implementation
% solution by x=b\A
y1=0;y2=0;y3=0;y4=0;y5=0;y6=0;y7=0;y8=0;
% solution by Gauss Seidel Method
% omitted
iters = 20000;
% solution by CBCD with block size 1
%[x1,y1] = CBCD_size1(A, b, d, lower, upper, iters);

%[x10,counter] = CoordDescentQPBox(A,b,lower, upper,lower);
% solution by RBCD with block size 1
%[x2,y2] = RBCD_size1(A, b, d, lower, upper, iters);
% solution by PDAL
%[x2, y2] = PDAL(A, b);
%[x02,y02,loop] = PDAL_mex(A,b,A'*A,A'*b,1/sqrt(max(eig(A'*A))));

% solution by CBCDmex with block size 1
 %[x3, y3] = CBCD_size3_gc(A, b, d, iters,1E-2,2,1,1);
[rx3, ry3] = RBCD_size3_gc(A, b, d, iters,1E-5,-1,1,0,1.0);
[rx2, ry2] = RBCD_size2_gc(A, b, d, iters,1E-5,-1,1,0,1.0);
[rx1, ry1] = RBCD_size1_gc(A, b, d, iters,1E-5,-1,1,0,1.0);
 %[x1, y1] = CBCD_size1_gc(A, b, d, iters,1E-2,2,1,1);
 %A = U'*A*U;A = sparse(A);
 %[x3, y3] = CBCD_size3_gc(A, b, d, iters,1E-2,2,1,1);
 [x2, y2] = CBCD_size1_gc(A, b, d, iters,1E-5,-1,1,0);
 %[x1, y1] = CBCD_size1_gc(A, b, d, iters,1E-2,2,1,1);
%[x3e, y3e] = CBCD_size1_mex_eigen(A, b, d, iters);
%x=CoordinateDescentQPBoxNonParallel(b,A);
% solution by RBCDmex with block size 1
%[x4, y4] = RBCD_size1_mex(A, b, d, lower, upper, iters);

% solution by CBCD with block size 2
%[x5, y5] = CBCD_size2(A, b, d, lower, upper, iters);
% solution by CBCDmex with block size 2
%[x70,y70] = CBCD_size2_9_mex_store(A, b, d, iters);
%[x7,y7] = CBCD_size2_9_mex(A, b, d, iters);
% solution by CBCD with block size 3
%[x8, y8] = CBCD_size3(A, b, d, lower, upper, iters);

%[x90,y90]  = CBCD_size3_mex_27_store(A, b, d, iters);
%[x9,y9]  = CBCD_size3_mex_27(A, b, d, iters);

%% sparse matrix implementation
% A = sparse(A);
%[x2ss, y2ss] = CBCD_size2_ss(A, b, d, iters);
%pause(1);
%[x3ss, y3ss] = CBCD_size3_ss(A, b, d, iters);
% [x1ss, y1ss] = CBCD_size1_mex_sparse(B, Bb, d, iters);
% pause(1);
% [x1ss, y1ss] = CBCD_size1_mex_sparse(C, Cb, d, iters);
% pause(1);
% [x1ss, y1ss] = CBCD_size1_mex_sparse(B, Bb, d, iters);




% [x1sp, y1sp] = CBCD_size1_mex_sparse(A, b, d, iters);
% [x2sp, y2sp] = CBCD_size2_9_mex_sparse(A, b, d, iters);
% [x3sp, y3sp] = CBCD_size3_mex_27_sparse(A, b, d, iters);

%% plot 
p = min([y1;y2;y3;y4;y5;y6;y7;y8]);
% figure(1),
% clf;
% semilogy(1:size(y1,1),y1-p,'LineWidth',2.5);
% hold on;
% semilogy(1:size(y2,1),y2-p,'r','LineWidth',2.5);
% hold on;
%semilogy(1:size(y3,1),y3-p,'g--','LineWidth',2.5);
% hold on;
% semilogy(1:size(y4,1),y4-p,'c','LineWidth',2.5);
%hold on;
% semilogy(1:size(y5,1),y5-p,'c--','LineWidth',2.5);
% hold on;
% semilogy(1:size(y6,1),y6-p,'m--','LineWidth',2.5);
% hold on;
% semilogy(1:size(y8,1),y8-p,'r--','LineWidth',2.5);
% hold on;
% legend
%legend('CBCD1','CBCD2','CBCD3');
%xlabel('#epoch','fontsize',16);ylabel('log10(f(x)-p*)','fontsize',16);
%set(gca,'fontsize',16);