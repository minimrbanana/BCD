%% compare the init value and the convergence speed
%% idx =200,300,400 for random matrix
addpath ../
exp_0=load('../result/EXP_idx_200/EXP.mat');
exp_5=load('../result/EXP_idx_300/EXP.mat');
exp_1=load('../result/EXP_idx_400/EXP.mat');
epoch0a=exp_0.EXP.epoch1;
epoch5a=exp_5.EXP.epoch1;
epoch1a=exp_1.EXP.epoch1;
epoch0c=exp_0.EXP.epoch3;
epoch5c=exp_5.EXP.epoch3;
epoch1c=exp_1.EXP.epoch3;
N=length(epoch0a);
kkt0a=zeros(N,1);
for i=1:N kkt0a(i)=exp_0.EXP.KKT_A{i}(1); end
kkt5a=zeros(N,1);
for i=1:N kkt5a(i)=exp_5.EXP.KKT_A{i}(1); end
kkt1a=zeros(N,1);
for i=1:N kkt1a(i)=exp_1.EXP.KKT_A{i}(1); end
kkt0c=zeros(N,1);
for i=1:N kkt0c(i)=exp_0.EXP.KKT_C{i}(1); end
kkt5c=zeros(N,1);
for i=1:N kkt5c(i)=exp_5.EXP.KKT_C{i}(1); end
kkt1c=zeros(N,1);
for i=1:N kkt1c(i)=exp_1.EXP.KKT_C{i}(1); end
figure(1),
d_kkta=kkt0a-kkt1a;
d_epocha=log(abs(epoch0a-epoch1a)).*sign(epoch0a-epoch1a);
plot(d_kkta,d_epocha,'.');
xlabel('\Deltakkt');
ylabel('\pmlog(\Deltaepoch)');
title('Initial KKT Condition and Convergence');
set(gca,'fontsize',18,'fontweight', 'bold');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_test_init/test_init_random_1.pdf');

product=d_kkta.*d_epocha;
p1=sum(product>0);
n1=sum(product<0);
figure(2),
d_kktc=kkt0c-kkt1c;
d_epochc=log(abs(epoch0c-epoch1c)).*sign(epoch0c-epoch1c);
plot(d_kktc,d_epochc,'.');
xlabel('\Deltakkt');
ylabel('\pmlog(\Deltaepoch)');
title('Initial KKT Condition and Convergence');
set(gca,'fontsize',18,'fontweight', 'bold');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_test_init/test_init_random_2.pdf');

product=d_kktc.*d_epochc;
p2=sum(product>0);
n2=sum(product<0);
%% compare init = 0.5
d1=epoch5a-epoch0a;sum(d1<0),d1(abs(d1)>500)=[];
d2=epoch5a-epoch1a;sum(d2<0),d2(abs(d2)>500)=[];
d3=epoch5c-epoch0c;sum(d3<0),d3(abs(d3)>500)=[];
d4=epoch5c-epoch1c;sum(d4<0),d4(abs(d4)>500)=[];
figure(12),title('Histogram of \Deltaepoch');
subplot(2,2,1),hist(d1,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
subplot(2,2,2),hist(d2,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
subplot(2,2,3),hist(d3,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
subplot(2,2,4),hist(d4,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_test_init/test_init_hist5randzoom.pdf');




%% idx =100,106,107 for 3-band matrix
exp_0=load('../result/EXP_idx_108/EXP.mat');
exp_5=load('../result/EXP_idx_107/EXP.mat');
exp_1=load('../result/EXP_idx_106/EXP.mat');
epoch0a=exp_0.EXP.epoch1;
epoch5a=exp_5.EXP.epoch1;
epoch1a=exp_1.EXP.epoch1;
epoch0c=exp_0.EXP.epoch3;
epoch5c=exp_5.EXP.epoch3;
epoch1c=exp_1.EXP.epoch3;
N=length(epoch0a);
kkt0a=zeros(N,1);
for i=1:N kkt0a(i)=exp_0.EXP.KKT_A{i}(1); end
kkt5a=zeros(N,1);
for i=1:N kkt5a(i)=exp_5.EXP.KKT_A{i}(1); end
kkt1a=zeros(N,1);
for i=1:N kkt1a(i)=exp_1.EXP.KKT_A{i}(1); end
kkt0c=zeros(N,1);
for i=1:N kkt0c(i)=exp_0.EXP.KKT_B{i}(1); end
kkt5c=zeros(N,1);
for i=1:N kkt5c(i)=exp_5.EXP.KKT_B{i}(1); end
kkt1c=zeros(N,1);
for i=1:N kkt1c(i)=exp_1.EXP.KKT_B{i}(1); end
figure(3),
d_kkta=kkt0a-kkt1a;
d_epocha=epoch0a-epoch1a;
plot(d_kkta,d_epocha,'.');
xlabel('\Deltakkt');
ylabel('\Deltaepoch');
title('Initial KKT Condition and Convergence');
set(gca,'fontsize',18,'fontweight', 'bold');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_test_init/test_init_band_1.pdf');

product=d_kkta.*d_epocha;
p3=sum(product>0);
n3=sum(product<0);
figure(4),
d_kktc=kkt0c-kkt1c;
d_epochc=epoch0c-epoch1c;
plot(d_kktc,d_epochc,'.');
xlabel('\Deltakkt');
ylabel('\Deltaepoch');
title('Initial KKT Condition and Convergence');
set(gca,'fontsize',18,'fontweight', 'bold');
grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_test_init/test_init_band_2.pdf');

product=d_kktc.*d_epochc;
p4=sum(product>0);
n4=sum(product<0);
%% compare init = 0.5
d5=epoch5a-epoch0a;sum(d5<0)
d6=epoch5a-epoch1a;sum(d6<0)
d7=epoch5c-epoch0c;sum(d7<0)
d8=epoch5c-epoch1c;sum(d8<0)
figure(6),title('Histogram of \Deltaepoch');
subplot(2,2,1),hist(d5,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
subplot(2,2,2),hist(d6,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
subplot(2,2,3),hist(d7,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
subplot(2,2,4),hist(d8,50);grid on;xlabel('\Deltaepoch');ylabel('Count');
set(gca,'fontsize',12,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_test_init/test_init_hist5.pdf');
