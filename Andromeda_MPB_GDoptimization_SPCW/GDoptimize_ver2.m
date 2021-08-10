%% 2017.10.10. Gradient descent optimization code // ver2
clear all
close all

Myname='Youn'
dates='20171014'
workname='SPCW_GDTest_test_v11'
%% Gradient Descent Optimzation by MPB simulation across the compute nodes in UQML Andromeda Cluster
currfolder=pwd;
fname='SPCW_ver2'
Myctlname=[fname '.ctl'];
Myoutname=[fname '.out'];
Myctldir=[currfolder '/SPCW'];
Andro_MPBdir='~/apps/mpb-1.5/bin/';
machinefiledir='machines';
workspace=['~/DATA/MPB/' dates '_' Myname '_' workname];

%% Constants and Unit conversion
c=299792458;
hbar=1.05457148e-34;
d2_freq=351.725; % THz
bluemagic_freq=c/(793.5*10^-9)*10^(-12);

%% Files preparation
if ~isdir(workspace)
    unix(['mkdir ' workspace])
end
unix(['cp ' Myctldir '/' Myctlname ' ' workspace '/' Myctlname]);
unix(['cp ' currfolder '/readBands.m ' workspace '/readBands.m']);
unix(['cp ' currfolder '/savetofile.m ' workspace '/savetofile.m']);
unix(['cp ' currfolder '/bandplot.m ' workspace '/bandplot.m']);
unix(['cp ' currfolder '/Modearea.m ' workspace '/Modearea.m']);
unix(['cp ' currfolder '/h5topng.m ' workspace '/h5topng.m']);
unix(['cp ' currfolder '/Fitband.m ' workspace '/Fitband.m']);
unix(['cp ' currfolder '/MPB_Run.m ' workspace '/MPB_Run.m']);
unix(['cp ' currfolder '/DATAProcess.m ' workspace '/DATAProcess.m']);

unix(['cp ' currfolder '/machine0 ' workspace '/machine0']);
unix(['cp ' currfolder '/machine1 ' workspace '/machine1']);
unix(['cp ' currfolder '/machine2 ' workspace '/machine2']);
unix(['cp ' currfolder '/machine3 ' workspace '/machine3']);
unix(['cp ' currfolder '/machine4 ' workspace '/machine4']);
unix(['cp ' currfolder '/machine5 ' workspace '/machine5']);
unix(['cp ' currfolder '/machine6 ' workspace '/machine6']);
cd(workspace)
%% Structure fixed parameters %%
thk=200;
wid=114;
r2=110;
r3=110;
s2=0;
FixedParam=[r2 thk wid r3 s2];

%% MPB AUX setting %%
nbands=24;
kstart=0.48;
kend=0.5;
kpts=10;
res=12;
Eta=[0.05 0.05 0.001 0.05 0.02];
Maxiter=40;
%% Initial Parameter set
a=374;
r1=110;
alpha=1.224;
offset=95;
s1=-5;
Paramset=zeros(5,Maxiter);
D=zeros(5,Maxiter);
Paramset(:,1)=Paramset(:,1)+[a;r1;alpha;offset;s1];
for i=1:Maxiter
    EpsParam{i}=[0.2 0.2 0.002 0.2 0.2].*exp(-i/6);
end
%% Run MPB simulation over the seven compute nodes
F=0.01;
Niter=1;
Nkpt=14;
Nband=20;
Nmachine=[0 1 2 3 4 5];
for p1234=1:Maxiter %F(Niter)>0.01 && Niter<Maxiter
    cd(workspace)
    CurrParam=Paramset(:,Niter);
%     if F<0.3
%         EpsParam=[1 1 0.05 1 1];
%     else
%         EpsParam=[2 2 0.05 2 2];
%     end
    RunParam{1}=CurrParam;
    Eps=EpsParam{Niter};
    Param=zeros(numel(Eps),1);
    for j=2:numel(Nmachine);        
        Param(j-1)=Paramset(j-1,Niter)+Eps(j-1);
        Param(1:end~=j-1)=Paramset(1:end~=j-1,Niter);
        RunParam{j}=Param;
    end
    %% N times MPB_Run to evalulate F(param(i)) & dF(param(i))
    parfor ii=1:numel(Nmachine)
        [data{ii}]=MPB_Run('mpbi-mpi',Myctlname,Myoutname,Andro_MPBdir,workspace,Nmachine(ii),FixedParam,RunParam{ii},nbands,kstart,kend,kpts,res,Niter);
    end
    FBand=vpa(data{1}.bands(end-Nkpt:end,Nband));
    DATAProcess(fname,FixedParam,RunParam{1},workspace,Nkpt,Niter)
    Flatness=0;
    a0=RunParam{1}(1)*10^(-9);
    SIunit_omega=(c/a0)*2*pi;  %%% frequency in THz
    trapomega=bluemagic_freq/(SIunit_omega./(2*pi).*10^(-12));
    d2omega=d2_freq/(SIunit_omega./(2*pi).*10^(-12));
    
    for j=1:length(FBand)-1
        Flatness=Flatness+vpa(abs(FBand(j)-FBand(j+1)));
    end
%     F(Niter+1)=10000*(8*Flatness+abs(vpa(FBand(end))-vpa(omega))+abs(vpa(data{1}.bands(end,22))-vpa(trapomega))).^2;
    F(Niter+1)=100*(100*Flatness);%+abs(vpa(FBand(end))-d2omega).^2);%+abs(vpa(data{1}.bands(end,22))-trapomega).^2);
    
    for jj=2:numel(Nmachine)
        DATAProcess(fname,FixedParam,RunParam{jj},workspace,Nkpt,Niter)
        Band{jj}=vpa(data{jj}.bands(end-Nkpt:end,Nband));
        DFlatness=zeros(5,1);
        ad(jj)=RunParam{jj}(1)*10^(-9);
        SIunit_omegad(jj)=(c/ad(jj))*2*pi;  %%% frequency in THz
        trapomegad(jj)=bluemagic_freq/(SIunit_omegad(jj)./(2*pi).*10^(-12));
        d2omegad(jj)=d2_freq/(SIunit_omegad(jj)./(2*pi).*10^(-12));
        
        for kk=1:length(Band{jj})-1
            DFlatness(jj-1)=DFlatness(jj-1)+vpa(abs(Band{jj}(kk)-Band{jj}(kk+1)));
        end
        DF(jj-1,Niter)=100*(100*DFlatness(jj-1));%+abs(Band{jj}(end)-d2omegad(jj)));%+abs(vpa(data{jj}.bands(end,22))-trapomegad(jj)));
        %             dF(jj,Niter)=DF(jj,Niter)*(vpa(DF(jj,Niter))-vpa(F(Niter+1)))/abs(vpa(DF(jj,Niter))-vpa(F(Niter+1)));
        dF(jj-1,Niter)=sqrt(F(Niter+1)).*(vpa(DF(jj-1,Niter))-vpa(sqrt(F(Niter+1))));%./(RunParam{jj}(jj-1)-RunParam{1}(jj-1));
    end
    for k=1:numel(Nmachine)-1
        Paramset(k,Niter+1)=Paramset(k,Niter)-vpa(dF(k,Niter))*Eta(k);
    end
%     fileID=fopen('F.txt','w');
%     fprintf(fileID,[num2str(Niter) ': ' num2str(F(Niter)) '\n ']);
    cd(workspace)
    H=figure(1579)
    h=plot(F,'ko-','linewidth',1.5)
    xlabel('Number of Iteration')
    ylabel('Objective function')
    saveas(H,'iterationcurve.png')
    save('Paramset.mat','Paramset')
    save('dF.mat','dF')
    save('F.mat','F')
    
    Niter=Niter+1;
end
cd(workspace)