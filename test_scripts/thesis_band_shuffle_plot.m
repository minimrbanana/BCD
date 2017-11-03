%% plot 3,5,7,9 shuffel band matrices
% parameters 
Fsize = 30; % font size
LineW = 6;  % line width


%% band 3
load /home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band3.mat

% #epochs to converge
figure(1);clf,
cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
plot(cover1,epoch1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD','Position',[0.5, 2.9E4]);
axis([0 1 0.5E4 2.6E4]);
%set(gcf,'Position',[232  246  560  440]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)*1.2];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/convergenceA3.pdf');

% spectral radius
figure(2);clf,
plot(cover1,Meig1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/spectralR3.pdf');

%% band 5
load /home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band5.mat

% #epochs to converge
figure(1);clf,
cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
plot(cover1,epoch1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD','Position',[0.5, 2.9E4]);
axis([0 1 0.5E4 2.6E4]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)*1.2];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/convergenceA5.pdf');

% spectral radius
figure(2);clf,
plot(cover1,Meig1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/spectralR5.pdf');

%% band 7
load /home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band7.mat

% #epochs to converge
figure(1);clf,
cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
plot(cover1,epoch1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD','Position',[0.5, 2.9E4]);
axis([0 1 0.5E4 2.6E4]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)*1.2];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/convergenceA7.pdf');

% spectral radius
figure(2);clf,
plot(cover1,Meig1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/spectralR7.pdf');

%% band 9
load /home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/band9.mat

% #epochs to converge
figure(1);clf,
cover1=0:0.1:1;cover2=0:0.1:1;cover3=0:0.1:1;
plot(cover1,epoch1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,epoch2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,epoch3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('#epoch to converge');
title('Convergence of CBCD','Position',[0.5, 2.9E4]);
axis([0 1 0.5E4 2.6E4]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)*1.2];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/convergenceA9.pdf');

% spectral radius
figure(2);clf,
plot(cover1,Meig1/Repeat,'ro-','LineWidth',LineW);hold on;
plot(cover2,Meig2/Repeat,'go-','LineWidth',LineW);hold on;
plot(cover3,Meig3/Repeat,'bo-','LineWidth',LineW);hold on;
legend('CBCD1','CBCD2','CBCD3','Location','southeast');
grid on;
xlabel('Permutation rate');
ylabel('max(|\lambda_{M}|)');
title('Spectral Radius of M');
axis([0 1 0.997 0.9993]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gca,'/home/yu/bcd/BCD/test_scripts/thesis_norm_of_M/spectralR9.pdf');

