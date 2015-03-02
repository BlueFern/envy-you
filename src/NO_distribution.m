% % Full activation:
% figure
% set(gcf, 'Position', [400 300 700 300]);
% y = [1.0615 0.2003; 0.6621 0.1438; 0.2627 0.08725; 0.2538 0.08598];
% y = [1.065 1.0615 0.2003; 0.6656 0.6621 0.1438; 0.2663 0.2627 0.08725; 0.2573 0.2538 0.08598];
% bar(y) %,0.4)
% Labels = {'NE', 'AC', 'SMC', 'EC'};
% set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
% ylabel('[NO] (\muM)')
% legend('maximum neuronal and endothelial NO release', 'neuronal', 'endothelial')


%% Physiological activation (resting state):
figure
set(gcf, 'Position', [400 300 700 300]);
% y = [1.0615 0.2003; 0.6621 0.1438; 0.2627 0.08725; 0.2538 0.08598];
y = [0.1814 0.1426 0.03814; 0.1249 0.08603 0.03814; 0.06841 0.02952 0.03815; 0.06714 0.02825 0.03815];
bar(y) %,0.4)
Labels = {'NE', 'AC', 'SMC', 'EC'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylabel('[NO] (\muM)')
legend('maximum neuronal and endothelial NO release', 'neuronal', 'endothelial')


%% During neuronal stimulation:
figure
set(gcf, 'Position', [400 300 700 300]);
y = [0.5151 0.4604 0.0508; 0.3326 0.2778 0.0508; 0.1501 0.0953 0.0508; 0.1460 0.0912 0.0508];
bar(y) %,0.4)
Labels = {'NE', 'AC', 'SMC', 'EC'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylabel('[NO] (\muM)')
legend('maximum neuronal and endothelial NO release', 'neuronal', 'endothelial')


%%
figure
set(gcf, 'Position', [400 300 700 300]);
Y = [0.1814 (1.6629-0.1814);
0.1426 (1.5076-0.1426);
0.0381 (0.1556-0.0381);

0.1249 (1.0651-0.1249);
0.0860 (0.9098-0.0860);
0.0381 (0.1556-0.0381);

0.0684 (0.4676-0.0684);
0.0295 (0.3121-0.0295);
0.0381 (0.1556-0.0381);

0.0671 (0.4542-0.0671);
0.0282 (0.2987-0.0282);
0.0381 (0.1556-0.0381)];

newY = reshape([reshape(Y,3,[]); zeros(1,numel(Y)/3)],[],2);  %# Add zeroes
bar(newY,'stacked');  %# Create a stacked histogram
shg
Labels = {'NE', 'AC', 'SMC', 'EC'};
set(gca,'XLim',[0 16],'XTick',2:4:14,'XTickLabel',Labels);  %# Modify axes

%% Scaled values:
% NOn = 1.0615;
% NOj = 0.2003; %0.08598;
% y = [1.0615/NOn 0.2003/NOj; 0.6621/NOn 0.1438/NOj; 0.2627/NOn 0.08725/NOj; 0.2538/NOn 0.08598/NOj];
