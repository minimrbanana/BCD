%% runtime of CBCD1 before and after reordering
% tested on four servers, average of 10 runs


%% BDC with block size 1 & 2 & 3

%% input
rng(1);

load /home/yu/datasets/SNAP/Social_networks/mat/facebook_combined.mat; iters=200000;
%load /home/yu/datasets/SNAP/Memetracker_and_Twitter/mat/higgs-social_network.edge.mat; iters=200;
%load /home/yu/datasets/SNAP/Collaboration_networks/mat/ca-CondMat.mat; iters=40000;
%load /home/yu/datasets/SNAP/Collaboration_networks/mat/ca-HepTh.mat; iters=200000;
%load /home/yu/datasets/SNAP/Location_social_networks/mat/loc-brightkite_edges.mat; iters=20000;
%load /home/yu/datasets/SNAP/Autonomous_systems_graphs/mat/as-caida20071105.mat; iters=40000;

d=size(A,1);
diag=sum(A);
A=A+speye(d)*1E-6;
xmin=rand(d,1)*2-0.5;
%xmin=rand(d,1)*0.8+0.1;
b=A*xmin;


%% true minimum and f min
% don't know
fmin=0.5*xmin'*A*xmin-b'*xmin;

init=0;

p1=symamd(A);
A1=A(p1,p1);b1=b(p1);
p2=symrcm(A);
A2=A(p2,p2);b2=b(p2);
t1=zeros(10,1);
t2=t1;t3=t1;
for i=1:10
    t0=tic;
    [cx1, cy1] = RBCD1(A, b, d, iters,1E-10,0,1,init);
    t1(i)=toc(t0);
    t0=tic;
    [cx2, cy2] = CBCD1(A1, b1, d, iters,1E-10,0,1,init);
    t2(i)=toc(t0);
    t0=tic;
    [cx3, cy3] = CBCD1(A2, b2, d, iters,1E-10,0,1,init);
    t3(i)=toc(t0);
    fprintf('*******\n i=%d\n*******\n',i);
end
[~, name] = system('hostname');
t_1=mean(t1);t_2=mean(t2);t_3=mean(t3);
st1=std(t1);st2=std(t2);st3=std(t3);
fid = fopen(['time_' name(1:3) '.txt'],'w');
fprintf(fid,'%.3f $\\pm$ %.3f\\\\ \n                %.3f $\\pm$ %.3f\\\\ \n                %.3f $\\pm$ %.3f} \n',...
    t_1,st1,t_2,st2,t_3,st3);
fclose(fid);
