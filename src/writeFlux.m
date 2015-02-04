function writeFlux(time,state)
% This function evaluates all the fluxxes and derivitives of the NVC Keet
% model and write them to a comma seperated file. This function can be
% implemented in the Outputfunction from the ODE solver where DEsyst is
% solved. This file uses flowrates and DEsyst to compute the actual flows
% and change of rates.
% The format of the CSVfile would be an array of the form:
% [AC(1:end), SMC(1:end), EC(1:end), dfdt(1:end), time(1);
%  AC(1:end), SMC(1:end), EC(1:end), dfdt(1:end), time(2);
%  ...]
numODE = length(state);
global csvfilename

    if size(time)> 0
        [NE,AC, SMC, EC] = all_fluxes(time(1),state);
        dfdt = DEsyst(time(1),state);
        dataNew = [NE,AC,SMC,EC,state((1:numODE),1)',dfdt((1:numODE),1)',time(1)...
            ,getRef(time(1),'ft')',getRef(time(1),'fluxft')',getRef(time(1),'rho')',getRef(time(1),'wss')',getRef(time(1),'Glu')'];

        try
            data = csvread(csvfilename);
            dataNew = [data;dataNew];
            csvwrite(csvfilename,dataNew);
        catch
            csvwrite(csvfilename,dataNew);
        end
    end
end