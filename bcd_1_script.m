%% BDC with block size 1
% Gauss Seidel Method

% input
% dimension and constraints
d = 100;
lower = zeros(d,1);
upper = ones(d,1);
% A and b
D = diag(rand(d,1));
U = orth(rand(d,d));
A = U' * D * U;
b = randn(d,1);

%solution
% solution by x=b\A
x = A\b;
% solution by Gauss Seidel Method
% omitted
iters = 99;
% solution by CBCD with block size 1
[x1,y1] = CBCD_size1(A, b, d, lower, upper, iters);

% solution by RBCD with block size 1
[x2,y2] = RBCD_size1(A, b, d, lower, upper, iters);

% solution by CBCDmex with block size 1
[x3, y3] = CBCD_size1_mex(A, b, d, lower, upper, iters);

% solution by CBCDmex with block size 1
[x4, y4] = RBCD_size1_mex(A, b, d, lower, upper, iters);

% plot 
p = min([y1;y2;y3;y4]);
figure(1),
clf;
semilogy(1:iters+1,y1-p);
hold on;
semilogy(1:iters+1,y2-p,'r');
hold on;
semilogy(1:iters+1,y3-p,'g--');
hold on;
semilogy(1:iters+1,y4-p,'c--');
%hold on;
% legend
legend('CBCD','RBCD','CBCD mex','RBCD mex');
xlabel('#epoch');ylabel('log10(f(x)-p*)');
% evaluation
diff1 = norm(x - x1,2);