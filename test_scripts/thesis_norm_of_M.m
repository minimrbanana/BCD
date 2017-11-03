%% load the norm and plot
% plot |M^{t}|_{2} with 0%,30%,60% 100% permutation
load Mnorm_0.mat
load Mnorm_30.mat
load Mnorm_60.mat
load Mnorm_100.mat


Fsize=20;


h41=figure(41);clf,
semilogy(1:size(M1norm_0,1),M1norm_0,'r','LineWidth',3);hold on;
semilogy(1:size(M2norm_0,1),M2norm_0,'g','LineWidth',3);hold on;
semilogy(1:size(M3norm_0,1),M3norm_0,'b','LineWidth',3);hold on;
grid on;
title('Matrix Norm of M');
xlabel('t');
ylabel('||M^{t}||_{2}');
axis([0 24000 1E-10 10]);
legend('CBCD1','CBCD2','CBCD3');
set(gca,'FontSize',Fsize,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/M1.pdf');

h42=figure(42);clf,
semilogy(1:size(M1norm_30,1),M1norm_30,'r','LineWidth',3);hold on;
semilogy(1:size(M2norm_30,1),M2norm_30,'g','LineWidth',3);hold on;
semilogy(1:size(M3norm_30,1),M3norm_30,'b','LineWidth',3);hold on;
grid on;
title('Matrix Norm of M');
xlabel('t');
ylabel('||M^{t}||_{2}');
axis([0 24000 1E-10 10]);
legend('CBCD1','CBCD2','CBCD3');
set(gca,'FontSize',Fsize,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/M2.pdf');

h43=figure(43);clf,
semilogy(1:size(M1norm_60,1),M1norm_60,'r','LineWidth',3);hold on;
semilogy(1:size(M2norm_60,1),M2norm_60,'g','LineWidth',3);hold on;
semilogy(1:size(M3norm_60,1),M3norm_60,'b','LineWidth',3);hold on;
grid on;
title('Matrix Norm of M');
xlabel('t');
ylabel('||M^{t}||_{2}');
axis([0 24000 1E-10 10]);
legend('CBCD1','CBCD2','CBCD3');
set(gca,'FontSize',Fsize,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/M3.pdf');

h44=figure(44);clf,
semilogy(1:size(M1norm_100,1),M1norm_100,'r','LineWidth',3);hold on;
semilogy(1:size(M2norm_100,1),M2norm_100,'g','LineWidth',3);hold on;
semilogy(1:size(M3norm_100,1),M3norm_100,'b','LineWidth',3);hold on;
grid on;
title('Matrix Norm of M');
xlabel('t');
ylabel('||M^{t}||_{2}');
axis([0 24000 1E-10 10]);
legend('CBCD1','CBCD2','CBCD3');
set(gca,'FontSize',Fsize,'fontweight', 'bold');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/M4.pdf');