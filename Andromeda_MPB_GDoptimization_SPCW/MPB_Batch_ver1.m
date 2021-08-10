%% 2017.10.02 Youn Lee
clear all
close all
Myname='Youn'
dates='20171024'
workname='SPCW_AlphaScan_Final'
%% MPB simulation across the compute nodes in UQML Andromeda Cluster
currfolder=pwd;
Myctlname='SPCW_ver2.ctl';
Myoutname='SPCW_ver2.out';
Myctldir=[currfolder '/SPCW'];
Andro_MPBdir='~/apps/mpb-1.5/bin/mpbi-mpi';
machinefiledir='machines';
workspace=['~/DATA/MPB/' dates '_' Myname '_' workname]

%% Constants and Unit conversion
c=299792458;
a0=373*10^-9;     % Lattice constant in m
ep0=8.85e-12;
hbar=1.05457148e-34;
h=6.626e-34;
MPB_unit=1;
d2_freq=351.725; % THz
bluemagic_freq=c/(793.5*10^-9)*10^(-12)
SIunit_omega=(c/a0)*2*pi*MPB_unit;  %%% frequency in THz
SIunit_energy=SIunit_omega*hbar;
SIunit_prop=(2*pi/a0);
SIunit_mass=hbar/(a0^2*2*pi*c/a0);
SIunit_k0=0.5*2*pi/a0;

%% Files preparation
if ~isdir(workspace)
    unix(['mkdir ' workspace])
end

unix(['cp ' Myctldir '/' Myctlname ' ' workspace '/' Myctlname])
unix(['cp ' currfolder '/readBands.m ' workspace '/readBands.m'])
unix(['cp ' currfolder '/savetofile.m ' workspace '/savetofile.m'])
unix(['cp ' currfolder '/bandplot.m ' workspace '/bandplot.m'])
unix(['cp ' currfolder '/Modearea.m ' workspace '/Modearea.m'])
unix(['cp ' currfolder '/h5topng.m ' workspace '/h5topng.m'])
unix(['cp ' currfolder '/Fitband.m ' workspace '/Fitband.m'])

unix(['cp ' currfolder '/machine0 ' workspace '/machine0'])
unix(['cp ' currfolder '/machine1 ' workspace '/machine1'])
unix(['cp ' currfolder '/machine2 ' workspace '/machine2'])
unix(['cp ' currfolder '/machine3 ' workspace '/machine3'])
unix(['cp ' currfolder '/machine4 ' workspace '/machine4'])
unix(['cp ' currfolder '/machine5 ' workspace '/machine5'])
unix(['cp ' currfolder '/machine6 ' workspace '/machine6'])
cd(workspace)
%% Structure parameters %%
a=373;
thk=200;
wid=114;
r1=110;
r2=110;
r3=110;
s1=0;
s2=0;
offset=100;
%alpha=1.25;
Nsq=4;

%% Parameter space set up
alpha=1.146:0.001:1.15;

%% MPB AUX setting %%
res=40;
kpts=100;
nbands=24;
kk=kpts+2;

%% Run MPB simulation over the seven compute nodes
Nmachine=[0 1 2 3 4];% 5];

parfor i=1:numel(Nmachine);
    Outputdir=[workspace '/res40_s1' num2str(s1(end)) '_r1' num2str(r1) '_r2' num2str(r2) ...
        '_s2' num2str(s2) '_alpha' num2str(alpha(i))];
    if ~isdir(Outputdir)
        unix(['mkdir ' Outputdir])
    end
    unix(['cp ' workspace '/' Myctlname ' ' Outputdir '/' Myctlname])
    unix(['cp ' workspace '/readBands.m ' Outputdir '/readBands.m'])
    unix(['cp ' workspace '/h5topng.m ' Outputdir '/h5topng.m'])
    unix(['cp ' workspace '/bandplot.m ' Outputdir '/bandplot.m'])
    unix(['cp ' workspace '/savetofile.m ' Outputdir '/savetofile.m'])
    unix(['cp ' workspace '/Fitband.m ' Outputdir '/Fitband.m'])
    unix(['cp ' workspace '/Modearea.m ' Outputdir '/Modearea.m'])
    unix(['cp ' workspace '/machine' num2str(Nmachine(i)) ' ' Outputdir '/machine' num2str(Nmachine(i))])
    
    cd(Outputdir)
    unix(['mpirun -np 12 --machinefile machine' num2str(Nmachine(i)) ' ' Andro_MPBdir ' ' ...
        'kpts=' num2str(kpts) ' nbands=' num2str(nbands) ' res=' num2str(res) ' s1=' num2str(s1) ...
        ' s2=' num2str(s2) ' a=' num2str(a) ' r1=' num2str(r1) ' r2=' num2str(r2) ' offset=' num2str(offset) ...
        ' alpha=' num2str(alpha(i)) ' ' Outputdir '/' Myctlname ' |tee ' Outputdir '/' Myoutname])
    
    %% Extract Data from OUT file and Save
    [data,bloc]=readBands(Myoutname);
    savetofile(data,'data.mat');
end

%% Data processing by Master node (bandplot, Intensity profile, Modearea, Delete H5 files)
for i=1:numel(Nmachine)
    Outputdir=[workspace '/res40_s1' num2str(s1(end)) '_r1' num2str(r1) '_r2' num2str(r2) ...
        '_s2' num2str(s2) '_alpha' num2str(alpha(i))];
    cd(Outputdir)
    %% Band plot
    bandplot(a*10^(-9));
    [Delta,Crrvature,Meff]=Fitband(20,10,a*10^(-9));
    %% Intensity profile
    %h5topng('SPCW_ver2',4,kk,20,20);
    %h5topng('SPCW_ver2',4,kk,22,22);
    %% Modearea
    [Aeff_even,Aeff_odd]=Modearea('SPCW_ver2',374e-9,20)
    Aeff20=Aeff_even;
    save('AeffProbe.mat','Aeff20')
    [Aeff_even,Aeff_odd]=Modearea('SPCW_ver2',374e-9,22)
    Aeff22=Aeff_odd;
    save('AeffTrap.mat','Aeff22')
    %% Delete H5file
    unix(['rm *.h5'])
    close all
end
