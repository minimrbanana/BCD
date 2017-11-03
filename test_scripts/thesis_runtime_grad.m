%% to test the run time of g=Ax before and after reordering and permutation
rng(1);
addpath ../
load /home/yu/datasets/SuiteSparse/circuit_simulation_problem/cvx/G2_circuit.mat; iters=500;
A=Problem.A;


%load /home/yu/datasets/SNAP/Networks_with_GT_communities/mat/com-amazon.ungraph.mat; iters=200;


d=size(A,1);
x=rand(d,1)*0.8+0.1;
b=A*x;

N=20;% repeat 50 times


p1=symamd(A);
p2=symrcm(A);
p3=randperm(d);
A1=A(p1,p1);b1=b(p1);
A2=A(p2,p2);b2=b(p2);
A3=A(p3,p3);b3=b(p3);

% plot structure of the 4 matrices
figure(1),set(gca,'fontsize',14);
subplot(1,4,1),spy(A);title('G2 circuit');
subplot(1,4,2),spy(A1);title('After AMD');
subplot(1,4,3),spy(A2);title('After RCM');
subplot(1,4,4),spy(A3);title('After Permutation');

% plot convergence
[~, cy1] = CBCD1(A, b, d, iters,1E-10,0,1,0.5);
[~, cy2] = CBCD1(A1, b1, d, iters,1E-10,0,1,0.5);
[~, cy3] = CBCD1(A2, b2, d, iters,1E-10,0,1,0.5);
[~, cy4] = CBCD1(A3, b3, d, iters,1E-10,0,1,0.5);

figure(2),set(gca,'fontsize',14);
subplot(2,2,1),
semilogy(0:iters-1,cy1,'r','linewidth',3);hold on;
semilogy(0:iters-1,cy2,'g','linewidth',2.5);hold on;
semilogy(0:iters-1,cy3,'b--','linewidth',2);hold on;
semilogy(0:iters-1,cy4,'c--','linewidth',1.5);hold on;
legend('A','A\_amd','A\_rcm','A\_perm');grid on;
xlabel('#epoch');ylabel('KKT Condition');
title('Convergence of Diverse Matrix A');

% repeat N times for runtime of solver
t1=zeros(N,1);t2=t1;t3=t1;t4=t1;
for i=1:N
    t0=tic;
    [~, cy1] = CBCD1(A, b, d, iters,1E-10,0,1,0.5);
    t1(i)=toc(t0);
    t0=tic;
    [~, cy2] = CBCD1(A1, b1, d, iters,1E-10,0,1,0.5);
    t2(i)=toc(t0);
    t0=tic;
    [~, cy3] = CBCD1(A2, b2, d, iters,1E-10,0,1,0.5);
    t3(i)=toc(t0);
    t0=tic;
    [~, cy4] = CBCD1(A3, b3, d, iters,1E-10,0,1,0.5);
    t4(i)=toc(t0);
end
subplot(2,2,2),
prc=prctile([t1,t2,t3,t4],[25 75],1);
mean_all = mean([t1,t2,t3,t4]);
errorbar(1,mean_all(1),prc(2,1)-mean_all(1),prc(1,1)-mean_all(1),'*');hold on;
errorbar(2,mean_all(2),prc(2,2)-mean_all(2),prc(1,2)-mean_all(2),'*');hold on;
errorbar(3,mean_all(3),prc(2,3)-mean_all(3),prc(1,3)-mean_all(3),'*');hold on;
errorbar(4,mean_all(4),prc(2,4)-mean_all(4),prc(1,4)-mean_all(4),'*');hold on;
legend('A','A\_amd','A\_rcm','A\_perm');grid on;
xlabel('Matrix Type');ylabel('Run Time(s)');
title('Run Time of Diverse Matrix A');


% repeat N times for runtime of A*x in matlab
t1=zeros(N,1);t2=t1;t3=t1;t4=t1;
for i=1:N
    t0=tic;
    y=A*x;
    t1(i)=toc(t0);
    t0=tic;
    y=A1*x;
    t2(i)=toc(t0);
    t0=tic;
    y=A2*x;
    t3(i)=toc(t0);
    t0=tic;
    y=A3*x;
    t4(i)=toc(t0);
end
subplot(2,2,3),
prc=prctile([t1,t2,t3,t4],[25 75],1);
mean_all = mean([t1,t2,t3,t4]);
errorbar(1,mean_all(1),prc(2,1)-mean_all(1),prc(1,1)-mean_all(1),'*');hold on;
errorbar(2,mean_all(2),prc(2,2)-mean_all(2),prc(1,2)-mean_all(2),'*');hold on;
errorbar(3,mean_all(3),prc(2,3)-mean_all(3),prc(1,3)-mean_all(3),'*');hold on;
errorbar(4,mean_all(4),prc(2,4)-mean_all(4),prc(1,4)-mean_all(4),'*');hold on;
legend('A','A\_amd','A\_rcm','A\_perm');grid on;
xlabel('Matrix Type');ylabel('Run Time(s)');
title('Run Time of A*x in MATLAB');


% repeat N times for runtime of A*x in C-mex
t1=zeros(N,1);t2=t1;t3=t1;t4=t1;
for i=1:N
    t0=tic;
    y=grad(A,x,d);
    t1(i)=toc(t0);
    t0=tic;
    y=grad(A1,x,d);
    t2(i)=toc(t0);
    t0=tic;
    y=grad(A2,x,d);
    t3(i)=toc(t0);
    t0=tic;
    y=grad(A3,x,d);
    t4(i)=toc(t0);
end
subplot(2,2,4),
prc=prctile([t1,t2,t3,t4],[25 75],1);
mean_all = mean([t1,t2,t3,t4]);
errorbar(1,mean_all(1),prc(2,1)-mean_all(1),prc(1,1)-mean_all(1),'*');hold on;
errorbar(2,mean_all(2),prc(2,2)-mean_all(2),prc(1,2)-mean_all(2),'*');hold on;
errorbar(3,mean_all(3),prc(2,3)-mean_all(3),prc(1,3)-mean_all(3),'*');hold on;
errorbar(4,mean_all(4),prc(2,4)-mean_all(4),prc(1,4)-mean_all(4),'*');hold on;
legend('A','A\_amd','A\_rcm','A\_perm');grid on;
xlabel('Matrix Type');ylabel('Run Time(s)');
title('Run Time of A*x in C-mex');

