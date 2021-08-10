%% batch 2017.10.30 
close all
clear;
% date='20171107'
% workname='SIL_x-pol_30degree_20nm_v4'
% workspace=['~/DATA/Lumerical/FDTD_DATA/' date '_' workname];
% if ~isdir(workspace)
%     mkdir(workspace)
% end
SimName='SideShiny';
%% Add Lumerical API to MATLAB
CurrentFolder=pwd;
path(path, 'C:\Program Files\Lumerical\FDTD\api\matlab')
%appclose(h);
h=appopen('fdtd');
%% Load Lumerical Header API
lheader(h,CurrentFolder)
%% Structure parameters
nSiO2= 1.45;
nSiN = 1.997;
%Time/Freq Domain info
fcs_d1=335.116e12;
fcs_d2=351.726e12;
a0=370e-9;
%% Dipole source parameters
wavelength=[600e-9 800e-9];
cc=299792458;
freqs=cc./wavelength;
cenfreq = (freqs(1)+freqs(2))/2;
bandwidth = freqs(1)-freqs(2);
freqhighres = 40; %nPoints over the frequency span

%% Simulation Parameters
meshacc = 4;
wSimVolx = 8.4e-6; %Simulation Volume : Squircle-W1 (80 unit cell)
wSimVoly = 4.2e-6;
wSimVolz = 3.2e-6;
dx=20e-9;
dy=dx;dz=dx;
timefactor = 10;

xSimMin=-wSimVolx;
xSimMax=wSimVolx;
ySimMin=-wSimVoly;
ySimMax=wSimVoly;
zSimMax=wSimVolz;
zSimMin=-wSimVolz;

%% Create the geometry
toggle=1;
% lgeometry(h,nSiN,dx,dy,dz,toggle)
STLfile=[CurrentFolder '/SPCW_v5_STL/SPCW_v5_40unitcell.stl']
lSTLimport(h,true,STLfile,nSiN,0,0,0,toggle)
%% Plane wave source
xpos_dip_up=0; ypos_dip_up=0; zpos_dip_up=3.0e-6;
xpos_dip_down=0; ypos_dip_down=0; zpos_dip_down=-3.0e-6;
Inject_thetad=[75,60,45,25,0,15,0];
Inject_phid=[90,0];
PolAngle=[90 0]; %90: y-pol, 0: x-pol
%PlaneWave(h,true,'PlaneWaveDOWN1',cenfreq,bandwidth,xpos_dip_down,ypos_dip_down,zpos_dip_down,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,'z',Inject_thetad(6),Inject_phid(1),PolAngle(1),'forward');
%PlaneWave(h,true,'PlaneWaveDOWN2',cenfreq+50e6,bandwidth,xpos_dip_down,ypos_dip_down,zpos_dip_down,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,'z',-Inject_thetad(6),Inject_phid(2),PolAngle(2),'forward');
PlaneWave(h,true,'PlaneWaveUP1',cenfreq+100e6,bandwidth,xpos_dip_up,ypos_dip_up,zpos_dip_up,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,'z',Inject_thetad(5),Inject_phid(2),PolAngle(2),'backward');
%PlaneWave(h,true,'PlaneWaveUP2',cenfreq+100e6,bandwidth,xpos_dip_up,ypos_dip_up,zpos_dip_up,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,'z',-Inject_thetad(2),Inject_phid(2),PolAngle(1),'backward');
%% Setting FDTD simulation area
SimX=0; SimY=0; SimZ=0;
AutoShutOff=1e-11;
Courant=1;
lFDTDcond(h,true,meshacc,timefactor,SimX,SimY,SimZ,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,dx,dy,dz,AutoShutOff,Courant);

%% Monitors
xpos=0;
ypos=0;
zpos=0;
lindex(h,true,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,0,0,0)
lFDTDmonitors(h,true,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,xpos,ypos,zpos,freqhighres)
%% Save and Run
lcommand(h,['save("' SimName '")'],true);
lcommand(h,'run',true);
lcommand(h,['save("' SimName '")'],true);
lcommand(h,'SILresult_X=getresult("SILDATA_X","E")',true);
lcommand(h,'SILresult_Y=getresult("SILDATA_Y","E")',true);
lcommand(h,'SILresult_Z=getresult("SILDATA_Z","E")',true);

lcommand(h,'n1 = getresult("index_monitor1","index")',true);
lcommand(h,'n2 = getresult("index_monitor2","index")',true);
lcommand(h,['matlabsave("nresult1",n1)'],true);
lcommand(h,['matlabsave("nresult2",n2)'],true);
lcommand(h,['matlabsave("SILresult_X",SILresult_X)'],true);
lcommand(h,['matlabsave("SILresult_Y",SILresult_Y)'],true);
lcommand(h,['matlabsave("SILresult_Z",SILresult_Z)'],true);
lcommand(h,'switchtolayout',true);
appclose(h)
%% Analysis
% unix(['mv ' CurrentFolder '/*.mat ' workspace ]) % move all MAT files to workspace
% unix(['mv ' CurrentFolder '/*.txt ' workspace ]) % move all MAT files to workspace
% unix(['mv ' CurrentFolder '/*.log ' workspace ]) % move all MAT files to workspace
% unix(['cp ' CurrentFolder '/Efieldplot_X_Normal.m ' workspace ]) % move Efieldplot.m file to workspace
% unix(['cp ' CurrentFolder '/Efieldplot_Y_Normal.m ' workspace ]) % move Efieldplot.m file to workspace
% unix(['cp ' CurrentFolder '/Efieldplot_Z_Normal.m ' workspace ]) % move Efieldplot.m file to workspace
% unix(['mv ' CurrentFolder '/' SimName '.fsp ' workspace]) % move fsp file to workspace
% cd(workspace)

nresult1=load(['nresult1.mat']);
nresult2=load(['nresult2.mat']);
nresult1=nresult1.n1;
nresult2=nresult2.n2;
SILresult_X=load(['SILresult_X.mat']);
SILresult_Y=load(['SILresult_Y.mat']);
SILresult_Z=load(['SILresult_Z.mat']);
SILresult_X=SILresult_X.SILresult_X;
SILresult_Y=SILresult_Y.SILresult_Y;
SILresult_Z=SILresult_Z.SILresult_Z;

lambdas=[613.5e-9,685.5e-9,793.5e-9];
ftrap=cc./lambdas;
for i=1:numel(ftrap)
    [Xplot,Yplot,Ex,Ey,Ez]=Efieldplot_X_Normal(SILresult_X,nresult1,ftrap(i));
    [Xplot,Yplot,Ex,Ey,Ez]=Efieldplot_Y_Normal(SILresult_Y,nresult1,ftrap(i));
    [Xplot,Yplot,Ex,Ey,Ez]=Efieldplot_Z_Normal(SILresult_Z,nresult2,ftrap(i));
end
