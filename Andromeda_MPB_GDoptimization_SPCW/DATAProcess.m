function DATAProcess(fname,FixedParam,CurrParam,workspace,Nkpt,Niter)
%% Data processing by Master node (bandplot, Intensity profile, Modearea, Delete H5 files)
Outputdir=[workspace '/Niter' num2str(Niter) '_s1' num2str(CurrParam(5),4) ...
    '_a' num2str(CurrParam(1),4) '_r1' num2str(CurrParam(2),4) '_alpha' num2str(CurrParam(3),4) ...
    '_offset' num2str(CurrParam(4),4)];
cd(Outputdir)
%% Band plot
bandplot(CurrParam(1)*10^(-9));
[Delta,Crrvature,Meff]=Fitband(20,Nkpt,CurrParam(1)*10^(-9));
%% Intensity profile
%h5topng('SPCW_ver2',4,kk,20,20);
%h5topng('SPCW_ver2',4,kk,22,22);
%% Modearea
% [Aeff_even,Aeff_odd]=Modearea(fname,CurrParam(1)*10^(-9),20)
% Aeff20=Aeff_even;
% save('AeffProbe.mat','Aeff20')
% [Aeff_even,Aeff_odd]=Modearea(fname,CurrParam(1)*10^(-9),22)
% Aeff22=Aeff_odd;
% save('AeffTrap.mat','Aeff22')
%% Delete H5file
unix(['rm *.h5'])
close all
