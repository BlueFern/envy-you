clear all
bet = 0.1300;
v_Ca3 = -27;
R_K = 12;

Ca = 0.2705; v = -35.66; 
for i = 1:20 cGMP2(i) = i; 
c_w(i)=i^10*0.0000000000006; y_Kactiv2(i) = ((Ca + c_w(i) )^2 / ( (Ca + c_w(i))^2 + bet*exp(-(v - v_Ca3)/R_K) ));
end
figure(13), plot(cGMP2,y_Kactiv2)
hold all

%%
clear all
bet = 0.1300;
v_Ca3 = -27;
R_K = 12;
c_w = 0; 
Ca = 0.2705;  
for c_w = [0,1]
    for i = -300:300 
        j = 301 +i;
        v(j) = i; 
        c_w ; 
        y_Kactiv2(j) =  ((Ca + c_w)^2 / ((Ca + c_w)^2 + bet*exp(-(v(j) - v_Ca3)/R_K) ));
    end
    figure(13), plot(v,y_Kactiv2)
    hold on
end

%%
clear all
bet = 0.1300;
v_Ca3 = -27;
R_K = 12;
c_w = 0; 
Ca = 0.2705;  
const1 =  0.0732;
for bet = [0.13,0.5]
    for i = -300:300 
        j = 301 +i;
        v(j) = i; 
        c_w ; 
        y_Kactiv2(j) = (const1) / ( (const1) + bet*exp(-(v(j) - v_Ca3)/R_K) );
    end
    figure(14), plot(v,y_Kactiv2)
    hold all
end


%% cGMP - c_w
%clear all

const1 = 0.0001; % shift in x-dir - the lower, the more right
const2 = 2; % steepness - the higher the steeper
const3 = 0; % shift in x-dir ??
const4 = 100000; % shift in x-dir - the higher, the more right

%best best
% const1 = 0.0001; % shift in x-dir - the lower, the more right
% const2 = 2; % steepness - the higher the steeper
% const3 = 0; % shift in x-dir ??
% const4 = 100000; % shift in x-dir - the higher, the more right

% best: 
% const1 = 0.0000001; % shift in x-dir - the lower, the more right
% const2 = 3; % steepness - the higher the steeper
% const3 = 1; % shift in x-dir ??
% const4 = 10000000; % shift in x-dir - the higher, the more right

% const1 = 0.0000000001; % shift in x-dir - the lower, the more right
% const2 = 5; % steepness - the higher the steeper
% const3 = 1; % shift in x-dir ??
% const4 = 100000000000; % shift in x-dir - the higher, the more right
% 
% const1 = 0.0000000000000001; % shift in x-dir - the lower, the more right
% const2 = 10; % steepness - the higher the steeper
% const3 = 1; % shift in x-dir ??
% const4 = 1000000000000000000000000; % shift in x-dir - the higher, the more right


for i = 0:20 
    j = 1 + i;
    cGMP(j) = i; 
    c_w(j) = const1 / (const1 + const4 * (exp(-cGMP(j) * const2 - const3)));
end
figure(14), plot(cGMP,c_w)
hold all
