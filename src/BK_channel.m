%% Dependencies of K_activation (=w_i)  on Ca_i and v_i

xlin = linspace(0.1,0.3,33);        % Ca_i - min / max
ylin = linspace(-60,-30,33);        % vi - min / max

c_w = 0;
R_K = 12;
bet = 0.1300;
v_Ca3 = -27;

[Cai,vi] = meshgrid(xlin,ylin);
Kact = ((Cai + c_w ).^2 ./ ( (Cai + c_w).^2 + bet.*exp(-(vi - v_Ca3)./R_K) ));
figure
surf(Cai,vi,Kact)
xlabel('Ca\_i')
ylabel('v\_i')
zlabel('K\_act') 
