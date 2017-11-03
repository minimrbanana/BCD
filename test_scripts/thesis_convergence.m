%% BDC with block size 1 & 2 & 3

%% input
%% BDC with block size 1 & 2 & 3

%% input
rng(1);
addpath ../
d = 30; 
e1 = ones(d,1);
lambda=1E-1;
lambdaI=speye(d)*lambda;
A=spdiags([-e1,2*e1,-e1],[-1,0,1],d,d);
A(1,d)=-1;A(d,1)=-1;

A=A+lambdaI;


% p1=symamd(A1);
% A2=A1(p1,p1);
% A0=A2;



xmin=zeros(d,1);
b0=xmin;

iters = 250;


%% plot 
figure(1),clf,
subplot(1,3,1),
imagesc(A),colorbar;
title('Matrix A');

[cx01, cy1] = CBCD1(A, b0, d, iters,1E-10,-2,1,-1);
[cx02, cy2] = CBCD2(A, b0, d, iters,1E-10,-2,1,-1);
[cx03, cy3] = CBCD3(A, b0, d, iters,1E-10,-2,1,-1);
subplot(1,3,2),
semilogy(1:size(cy1,1),cy1,'r','LineWidth',2.5);
hold on;
semilogy(1:size(cy2,1),cy2,'g','LineWidth',2.5);
hold on;
semilogy(1:size(cy3,1),cy3,'b','LineWidth',2.5);
hold on;
legend('CBCD1','CBCD2','CBCD3');
grid on;
title('Convergence');
xlabel('t');
ylabel('KKT Condition');
axis([0 250 1E-11 1]);
%% analysis
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

block1=diag(ones(d,1));

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
I=eye(d);
AL1=A.*L1;
M1=I-AL1\A;
AL2=A.*L2;
M2=I-AL2\A;
AL3=A.*L3;
M3=I-AL3\A;


Norm1=zeros(iters,1);
Norm2=Norm1;
Norm3=Norm1;
temp1=M1;
temp2=M2;
temp3=M3;
for i=1:iters
    Norm1(i)=norm(temp1);
    Norm2(i)=norm(temp2);
    Norm3(i)=norm(temp3);
    temp1=temp1*M1;
    temp2=temp2*M2;
    temp3=temp3*M3;
end


subplot(1,3,3),
semilogy(1:iters,Norm1,'r','LineWidth',2.5);hold on;
semilogy(1:iters,Norm2,'g','LineWidth',2.5);hold on;
semilogy(1:iters,Norm3,'b','LineWidth',2.5);hold on;
legend('CBCD1','CBCD2','CBCD3');
grid on;
title('Matrix Norm of M');
xlabel('t');
ylabel('||M^{t}||');
axis([0 250 1E-30 1]);



figure(2),clf,

subplot(1,3,1),
spy(M1);
title('Structure of M_1');

subplot(1,3,2),
spy(M2);
title('Structure of M_2');

subplot(1,3,3),
spy(M3);
title('Structure of M_3');


eig(M1)