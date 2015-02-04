clear all
close all
%% Oxygen:
K_mO2_j = 7.7; % (11.25) uM Chen2006
K_mO2_n = 243; % 140-145 uM Chen2007

for i = 1:200
    oxygen(i) = i;
    P_NO_j1(i) = 0.055 * oxygen(i) / (K_mO2_j + oxygen(i));
    P_NO_n1(i) = 3.8 * oxygen(i) / (K_mO2_n + oxygen(i));
end

figure; plot(oxygen,P_NO_j1)
xlabel('[O_2] in EC [uM]')
ylabel('NO production rate in EC [uM/s]')
figure; plot(oxygen,P_NO_n1)
xlabel('[O_2] in NC [uM]')
ylabel('NO production rate in NE [uM/s]')

%% L-Arginine:

K_mArg_j = 1.9; % (0.9-2.9, Chen2006) 
K_mArg_n = 2.3; % Babu1999 (at pH 7.4)

for i = 1:100
    arginine(i) = i;
    P_NO_j2(i) = 0.055 * arginine(i) / (K_mArg_j + arginine(i));
    P_NO_n2(i) = 3.8 * arginine(i) / (K_mArg_n + arginine(i));
end

figure; plot(arginine,P_NO_j2)
xlabel('[Arg] in EC [uM]')
ylabel('NO production rate in EC [uM/s]')
figure; plot(arginine,P_NO_n2)
xlabel('[Arg] in NC [uM]')
ylabel('NO production rate in NE [uM/s]')

%% both:

[x_arg,y_oxy] = meshgrid(0:3:100,0:3:200);
z_PNO = 0.055 .* x_arg ./ (K_mArg_j + x_arg) .* y_oxy ./ (K_mO2_j + y_oxy);
figure
surface(x_arg,y_oxy,z_PNO)
view(3)
xlabel('[Arg]')
ylabel('[O_2]')
zlabel('P_{NO} endothelial cell')

z_PNO = 3.8 .* x_arg ./ (K_mArg_n + x_arg) .* y_oxy ./ (K_mO2_n + y_oxy);
figure
surface(x_arg,y_oxy,z_PNO)
view(3)
xlabel('[Arg]')
ylabel('[O_2]')
zlabel('P_{NO} neuron')