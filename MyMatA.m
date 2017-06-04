function A = MyMatA(eidx, d)

lambda=1E-6;

switch eidx
    case 1
        error('matrix A is defined!');
    case 2
        error('matrix A is defined!');
    case 3
        % A has block size2 on the diagonal
        % with noise off the diagonal
        e0 = [1;0];
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,e1*0,-[0;e1(1:end-1)]],[-1,0,1],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 4
        % A has block size2 on the diagonal with shift
        % with noise off the diagonal
        e0 = [0;1]; % change from [0,1] to [1,0] as shift
        e1 = e0(:,ones(ceil(d/2),1));
        e1 = reshape(e1 ,numel(e1),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,e1*0,-[0;e1(1:end-1)]],[-1,0,1],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 5
        error('matrix A is defined!');
    case 6
        error('matrix A is defined!');
    case 7
        % A has block size3 on the diagonal
        % with noise off the diagonal
        e0 = [1;1;0];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;0;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e2,-e1,e1*0,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],...
            [-2,-1,0,1,2],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 8
        % A has block size3 on the diagonal with shift
        % with noise off the diagonal
        e0 = [0;1;1];
        e1 = e0(:,ones(ceil(d/3),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [0;1;0];
        e2 = e0(:,ones(ceil(d/3),1));
        e2 = reshape(e2 ,numel(e2),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e2,-e1,e1*0,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)]],...
            [-2,-1,0,1,2],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 9
        error('matrix A is defined!');
    case 10
        error('matrix A is defined!');        
    case 11
        % A has block size4 on the diagonal WITH noise
        e0 = [1;1;1;0];
        e1 = e0(:,ones(ceil(d/4),1));
        e1 = reshape(e1 ,numel(e1),1);
        e0 = [1;1;0;0];
        e2 = e0(:,ones(ceil(d/4),1));
        e2 = reshape(e2 ,numel(e2),1);
        e0 = [1;0;0;0];
        e3 = e0(:,ones(ceil(d/4),1));
        e3 = reshape(e3 ,numel(e3),1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e3,-e2,-e1,-[0;e1(1:end-1)],-[0;0;e2(1:end-2)],...
            -[0;0;0;e3(1:end-3)]],[-3,-2,-1,1,2,3],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 100
        error('matrix A is defined!');
    case 101
        error('matrix A is defined!');
    case 102
        error('matrix A is defined!');
    case 103
        % A is a tri-diagonal matrix WITH noise
        e1 = ones(d,1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,e1*0,-e1],[-1,0,1],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 104
        % A is a 5-band matrix WITH noise
        e1 = ones(d,1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,-e1,-e1,-e1],[-2,-1,1,2],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A); 
    case 105
        % A is a 7-band matrix WITH noise
        e1 = ones(d,1);
        A = sprandsym(d,3/d,0.5,1);
        A_off = -A./(A+eps);
        A = spdiags([-e1,-e1,-e1,-e1,-e1,-e1],[-3,-2,-1,1,2,3],A_off);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);          
    case 200
        % A is a random symmetric matrix 
        % sparsity 3/d
        A = sprandsym(d,3/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 201
        % A is a random symmetric matrix 
        % sparsity 5/d
        A = sprandsym(d,5/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 202
        % A is a random symmetric matrix 
        % sparsity 7/d
        A = sprandsym(d,7/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 203
        % A is a random symmetric matrix 
        % sparsity 9/d
        A = sprandsym(d,9/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 204
        % A is a random symmetric matrix 
        % sparsity 11/d
        A = sprandsym(d,11/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 205
        % A is a random symmetric matrix 
        % sparsity 13/d
        A = sprandsym(d,13/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 206
        % A is a random symmetric matrix 
        % sparsity 15/d
        A = sprandsym(d,15/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 207
        % A is a random symmetric matrix 
        % sparsity 17/d
        A = sprandsym(d,17/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 208
        % A is a random symmetric matrix 
        % sparsity 30/d
        A = sprandsym(d,30/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 209
        % A is a random symmetric matrix 
        % sparsity 40/d
        A = sprandsym(d,40/d,0.5,1);
        A = spdiags(zeros(d,1),0,A);
        A = -A./(A+eps);
        diagonal = -sum(A)+lambda;
        diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
        A = spdiags(diagonal',0,A);
    case 900
        error('matrix A is defined!');
    case 901
        error('matrix A is defined!');
    case 902
        error('matrix A is defined!');
    case 903
        c = 0.8;
        A = ones(d,d)*c;
        A = A + diag(ones(d,1)*(1-c));
        [U,~,~] = svd(1-2*rand(d,d));
        A = U'*A*U;
        A = sparse(A);
    
end


end