value=0;
while value==0
    save1 = input('Would you like to save all figures and parameters in a separate folder? (yes [y] / no [n]): ', 's');
    if strcmp(save1,'y') || strcmp(save1,'yes')
        fprintf('Please wait!');
        cl = clock;
        filename = ['Results-',date,'-',num2str(cl(4),'%.2d'),num2str(cl(5),'%.2d')];
        mkdir(filename);
        cd(filename);
        saveas(figure(1),'1 Input_Signal.fig');  %or fullfile()
        saveas(figure(2),'2 SMC_Fluxes.fig');
        saveas(figure(3),'3 EC_Fluxes.fig');
        saveas(figure(4),'4 DFDT.fig');
        saveas(figure(5),'5 Myosin-Crossbridge_Model.fig');
        saveas(figure(6),'6 Neuron_to_Radius.fig');
        saveas(figure(7),'7 Calcium_and_Radius.fig');
        
        file1=fopen('Parameters.txt','w');
        fprintf(file1, ['Parameters of "' filename '":\r\n\r\n']);
        fprintf(file1, 't_start      = %.0f      (s)\r\n', t_start);
        fprintf(file1, 't_end        = %.0f    (s)\r\n', t_end);
        fprintf(file1, 'startpulse   = %.0f    (s)\r\n', startpulse);
        fprintf(file1, 'lengthpulse  = %.0f    (s)\r\n', lengthpulse);
        fprintf(file1, 'CASE         = %.0f      (dim.less)\r\n', CASE);
        fprintf(file1, 'J_PLC        = %.3f  (muM s-1)\r\n', J_PLC);
        fprintf(file1, 'C_Hillmann   = %.0f      (dim.less)\r\n', C_Hillmann);
        fprintf(file1, 'stretch_ch   = %3s    (dim.less)\r\n', stretch_ch);
        fprintf(file1, 'only_Koenig  = %3s    (dim.less)\r\n', only_Koenig);
        
        
        fclose(file1);
        clc; fprintf('Figures and parameters have been saved in the folder "%s".\n', filename);
        value=1;
        cd ..
    elseif strcmp(save1,'n') || strcmp(save1,'no')
        clc; fprintf('Figures and parameters have not been saved.\n');
        value=1;
    else
        clc; fprintf('Warning: "%s" is not a valid input!\n', save1);
    end
end