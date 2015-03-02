
for i = 1:200
    j = i-100;
    x_vi(i) = j;
    y_VOOC(i) = G_Ca * ( x_vi(i) - v_Ca1) / ( 1 + exp( - ( x_vi(i) - v_Ca2 ) / ( R_Ca ) ) );
end
figure, plot(x_vi,y_VOOC)
    


%%
