function [B,R] = RCM(A)
%% Reverse CutHill-McKee Algorithm
% input matrix
%A = [1 0 0 0 1 0 0 0 ;0 1 1 0 0 1 0 1;0 1 1 0 1 0 0 0 ;0 0 0 1 0 0 1 0;1 0 1 0 1 0 0 0;0 1 0 0 0 1 0 1 ;0 0 0 1 0 0 1 0 ;0 1 0 0 0 1 0 1];
assert(size(A,1)==size(A,2));
assert(norm(A-A','fro')<1e-14);
disp('start RCM');
% dim of input
d = size(A,1);
% init output index
R = zeros(d,1);
S = sum(A);
[D, index]=sort(S,'ascend');
index_left = index;% this is for start new dis-connected graphs
% first element in R and Q
R(1) = index(1);
index_left(1)=[];
nextp = 2; % next position in R
Q = find(A(:,R(1))~=0);
Q = sort(Q,'ascend');
while nnz(R)<d % loop until R is filled with D numbers
    while ~isempty(Q) & find(R==Q(1))>0
        Q = Q(2:end);
    end
    if ~isempty(Q)
        % put first element in Q to R and remove it in Q
        R(nextp) = Q(1);
        %R'
        index_left(index_left==Q(1))=[];
        nextp=nextp+1;
        Qplus = find(A(:,Q(1))~=0);
        Q = [Q(2:end);Qplus];
    end
    % if Q is empty, i.e. dis-connected graph exists, pick new start
    if isempty(Q)
        R(nextp) = index_left(1);
        %R'
        Q = find(A(:,index_left(1))~=0);
        index_left(1)=[];
        nextp=nextp+1;
    end
    if mod(nextp,1000)==0
        disp('...');
    end
end
% reverse R
R = flipud(R);
B = A(R,:);
B = B(:,R);
disp('end RCM');
end



