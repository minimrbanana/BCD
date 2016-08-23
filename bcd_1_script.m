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
iters = 500;
x1 = GaussSeidel(A, b, d, lower, upper, iters);
% solution by CBCD with block size 1
iters = 500;
x2 = CBCD_size1(A, b, d, lower, upper, iters);
% solution by RBCD with block size 1
iters = 500;
x3 = RBCD_size1(A, b, d, lower, upper, iters);
% legend
legend('G-S','CBCD','RBCD');
% evaluation
diff1 = norm(x - x1,2);
diff2 = norm(x - x2,2);