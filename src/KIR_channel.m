%% Dependencies of J_KIR_i on K_p and v_i
xlin = linspace(0,15,33);           % Kp - min / max
ylin = linspace(-80,-20,33);  % vi - min / max

[Kp,vi] = meshgrid(xlin,ylin);
vKIR = z_1 .* K_p./unitcon + z_2;
GKIR = exp(z_5 .* vi + z_3 .* Kp./unitcon + z_4 );        
JKIR = F_il .* GKIR .*(vi-vKIR)./gam;    

figure
surf(Kp,vi,JKIR)
xlabel('K\_p')
ylabel('v\_i')
zlabel('J\_KIR')
% 
% vKIR = z_1 .* K_p./unitcon + z_2;
% GKIR = exp(z_5 .* vi + z_3 .* Kp./unitcon + z_4 );        
% JKIR = F_il./gam .* GKIR .*(vi-vKIR);    