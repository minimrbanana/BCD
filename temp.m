% get the matrix
d = 5000;
RC=[ones(d/4,1);ones(d/4,1)*10;ones(d/4,1)*20;ones(d/4,1)*30];
A = sprandsym(d,5/d,RC);
iters=2000000;

% set outputs
n=10;
epochs2=zeros(1,n);
epochs2h=zeros(1,n);
epochs2u=zeros(1,n);
epochs3=zeros(1,n);
epochs3h=zeros(1,n);
epochs3u=zeros(1,n);
T=zeros(4,n);
for loop=1:n
    b = randn(d,1);
    profile on;
    [~,r3h]=RBCD_size3_gc_h(A, b, d, iters,1E-5,0,1,0,1.0);
    [~,r3u]=RBCD_size3_gc_u(A, b, d, iters,1E-5,0,1,0,1.0);
    [~,r3 ]=RBCD_size3_gc(A, b, d, iters,1E-5,0,1,0,1.0);
    [~,r2h]=RBCD_size2_gc_h(A, b, d, iters,1E-5,0,1,0,1.0);
    [~,r2u]=RBCD_size2_gc_u(A, b, d, iters,1E-5,0,1,0,1.0);
    [~,r2 ]=RBCD_size2_gc(A, b, d, iters,1E-5,0,1,0,1.0);
    p=profile('info');
    profile off;
    % time 
    for i=1:size(p.FunctionTable,1)
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size3_gc_h')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(1,loop)= p.FunctionTable(i,1).TotalTime;
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size3_gc_u')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(2,loop)= p.FunctionTable(i,1).TotalTime;
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size3_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(3,loop)= p.FunctionTable(i,1).TotalTime;
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size2_gc_h')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(4,loop)= p.FunctionTable(i,1).TotalTime;
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size2_gc_u')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(5,loop)= p.FunctionTable(i,1).TotalTime;
        end
        if strcmp(p.FunctionTable(i,1).FunctionName , 'RBCD_size2_gc')
            fprintf('Function: %s\n',p.FunctionTable(i,1).FunctionName);
            T(6,loop)= p.FunctionTable(i,1).TotalTime;
        end
    end
    epochs2(loop)=size(r2,1);
    epochs2h(loop)=size(r2h,1);
    epochs2u(loop)=size(r2u,1);
    epochs3(loop)=size(r3,1);
    epochs3h(loop)=size(r3h,1);
    epochs3u(loop)=size(r3u,1);
end
h2=mean(epochs2h./epochs2);
h3=mean(epochs3h./epochs3);
u2=mean(epochs2u./epochs2);
u3=mean(epochs3u./epochs3);
t3h=mean(T(1,:)./T(3,:));
t3u=mean(T(2,:)./T(3,:));
t2h=mean(T(4,:)./T(6,:));
t2u=mean(T(5,:)./T(6,:));
fprintf('h2=%.5f, h3=%.5f\n',h2,h3);
fprintf('u2=%.5f, u3=%.5f\n',u2,u3);
fprintf('t3h=%.5f, t3u=%.5f\n',t3h,t3u);
fprintf('t2h=%.5f, t2u=%.5f\n',t2h,t2u);






