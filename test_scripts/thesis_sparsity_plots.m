%% some plots
Fsize = 18;
% plot the speed up with increasing of sparsity
res2o1 = zeros(1,9);
res3o1 = zeros(1,9);
i=200;
filename='../result/EXP_idx_300/EXP.mat';
load(filename);
res2o1(i-199) = EXP.C_ep2_OVER_ep1;
res3o1(i-199) = EXP.C_ep3_OVER_ep1;

for i=201:208
    filename=['../result/EXP_idx_' num2str(i) '/EXP.mat'];
    load(filename);
    res2o1(i-199) = EXP.C_ep2_OVER_ep1;
    res3o1(i-199) = EXP.C_ep3_OVER_ep1;
end
sparsity = [3,5,7,9,11,13,15,17,30]/3000;
%res2o1 = [0.8519,0.8716,0.7593,0.9052,0.8755,0.8222,0.8816,0.8682,0.8736,0.9969];
%res3o1 = [0.7536,0.6469,0.7391,0.8523,0.8633,0.8960,0.8421,0.7767,0.8761,0.8941];
figure(1),
plot(sparsity, res2o1,'r','LineWidth',2);hold on; 
plot(sparsity, res3o1,'b','LineWidth',2);
xlabel('Sparsity: s/dim');
ylabel('Speed up');
legend('#CBCD2/#CBCD1','#CBCD3/#CBCD1','Location','SouthEast');
axis([2/3000,30/3000,0.7 1]);grid on;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'FontSize',Fsize,'fontweight', 'bold');
saveas(gcf,'sparsity.pdf');
