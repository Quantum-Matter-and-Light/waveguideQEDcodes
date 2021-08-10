%% 2017.10.08. Gradient descent optimization code
clear all
close all

Myname='Youn'
dates='20171010'
workname='SPCW_GradientDescentOpt'
%% Gradient Descent Optimzation by MPB simulation across the compute nodes in UQML Andromeda Cluster
currfolder=pwd;
fname='SPCW_ver3'
Myctlname=[fname '.ctl'];
Myoutname=[fname '.out'];
Myctldir=[currfolder '/SPCW'];
Andro_MPBdir='~/apps/mpb-1.5/bin/';
machinefiledir='machines';
workspace=['~/DATA/MPB/' dates '_' Myname '_' workname];

%% Constants and Unit conversion
c=299792458;
a0=370*10^-9;     % Lattice constant in m
hbar=1.05457148e-34;
d2_freq=351.725; % THz
bluemagic_freq=c/(793.5*10^-9)*10^(-12);
SIunit_omega=(c/a0)*2*pi;  %%% frequency in THz
SIunit_energy=SIunit_omega*hbar;
SIunit_prop=(2*pi/a0);
SIunit_mass=hbar/(a0^2*2*pi*c/a0);
SIunit_k0=0.5*2*pi/a0;
trapomega=bluemagic_freq/(SIunit_omega./(2*pi).*10^(-12));
omega=d2_freq/(SIunit_omega./(2*pi).*10^(-12));
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
FixedParam=[r2 thk wid r3];

%% MPB AUX setting %%
nbands=24;
kstart=0.42;
kend=0.5;
kpts=16;
res=16;
Eta=[20 20 0.6 20 20 20];
Maxiter=30;
%% Initial Parameter set
a=372;
r1=110;
alpha=1.2;
offset=130;
s1=-10;
s2=-10;
Paramset=zeros(6,Maxiter);
D=zeros(6,Maxiter);
Paramset(:,1)=Paramset(:,1)+[a;r1;alpha;offset;s1;s2];
%% Run MPB simulation over the seven compute nodes
F=1;
Niter=1;
Nkpt=10;
Nband=20;
Nmachine=[0 1 2 3 4 5];
while F(Niter)>0.1 && Niter<Maxiter
    %% MPB_Run 
    CurrParam=Paramset(:,Niter);
    [Data]=MPB_Run(Myctlname,Myoutname,Andro_MPBdir,workspace,Nmachine(1),FixedParam,CurrParam,nbands,kstart,kend,kpts,res,Niter);
    DATAProcess(fname,FixedParam,CurrParam,workspace,Nkpt,Niter);
    %% Evaluates Objective function "F(paramset)"
    FBand=Data.bands(end-Nkpt:end,Nband);
    Flatness=0;
    for j=1:length(FBand)-1
        Flatness=Flatness+abs(vpa(FBand(j))-vpa(FBand(j+1)));
    end
    F(Niter+1)=100*(10*Flatness+abs(vpa(FBand(end))-vpa(omega))+abs(vpa(Data.bands(end,22))-vpa(trapomega)));
    if F(Niter+1)>0.1
        cd(workspace)
%        Niter=Niter+1;
        %% N times MPB_Run to evalulate dF(param(i))
        DParam=zeros(numel(Nmachine),1);
        if F<0.3
           EpsParam=[1 1 0.02 1 1 1];
        else
           EpsParam=[2 2 0.05 2 2 2];
        end
        for i=1:numel(Nmachine);
            DParam(i)=Paramset(i,Niter)+EpsParam(i);
            DParam(1:end~=i)=Paramset(1:end~=i,Niter);
            DeltaParam{i}=DParam;
        end        
        parfor ii=1:numel(Nmachine)
            [data{ii}]=MPB_Run(Myctlname,Myoutname,Andro_MPBdir,workspace,Nmachine(ii),FixedParam,DeltaParam{ii},nbands,kstart,kend,kpts,res,Niter);
        end
        for jj=1:numel(Nmachine)        
            DATAProcess(fname,FixedParam,DeltaParam{jj},workspace,Nkpt,Niter)
            Band=data{jj}.bands(end-Nkpt:end,Nband);
            DFlatness=0;
            for kk=1:length(Band)-1
                DFlatness=DFlatness+abs(vpa(Band(i))-vpa(Band(i+1)));
            end
            DF(jj,Niter)=100*(10*Flatness+abs(Band(end)-omega)+abs(vpa(data{jj}.bands(end,22))-vpa(trapomega)));
%             dF(jj,Niter)=DF(jj,Niter)*(vpa(DF(jj,Niter))-vpa(F(Niter+1)))/abs(vpa(DF(jj,Niter))-vpa(F(Niter+1)));
            dF(jj,Niter)=vpa(DF(jj,Niter))-vpa(F(Niter+1));
            Paramset(jj,Niter+1)=Paramset(jj,Niter)-vpa(dF(jj,Niter))*Eta(jj);
        end
    end
    Niter=Niter+1;
    cd(workspace)
end
