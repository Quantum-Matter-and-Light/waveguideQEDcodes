%% batch 2017.11.01 revised by Youn
close all
clear;
date='20171103'
workname='GF_60unitcell'
workspace=['E:/Photonic Crystal Codes/Lumerical Green Function Code/Youn/' date '_' workname];
if ~isdir(workspace)
    mkdir(workspace)
end
SimName='GreenFunction';
for id=1:1
    if id==2
        dipoleaxis='x'
    elseif id==1
        dipoleaxis='y'
    elseif id==3
        dipoleaxis='z'
    end
    cc=299792458; %speed of light
    %% Add Lumerical API to MATLAB
    CurrentFolder=pwd;
    path(path, 'C:\Program Files\Lumerical\FDTD\api\matlab')
    %appclose(h);
    h=appopen('fdtd');
    %% Load Lumerical Header API
    lheader(h,CurrentFolder)
    %% Structure parameters
    nSiN = 1.997;
    fcs_d1=335.116e12;
    fcs_d2=351.726e12;
    a0=370e-9;
    %% Dipole source parameters
    freq = 350e12;;
    bandwidth = 10e12;
    freqhighres = 10000; %nPoints over the frequency span
    %% Simulation Parameters
    meshacc = 4;
    wSimVolx = 1.8*8.8e-6;
    wSimVoly = 4.4e-6;
    wSimVolz = 0.8e-6;
    dx=40e-9;
    dy=dx;dz=dx;
    timefactor =30000;
    
    xSimMin=-wSimVolx;
    xSimMax=wSimVolx;
    ySimMin=-wSimVoly;
    ySimMax=wSimVoly;
    zSimMax=wSimVolz;
    zSimMin=-wSimVolz;    
    %% Create the geometry
    toggle=1;
%     lgeometry(h,nSiN,dipoleaxis,dx,dy,dz,toggle)    
    xpos_SPCW=0; ypos_SPCW=0; zpos_SPCW=0;
    xpos_dip=0; ypos_dip=0; zpos_dip=0;
    SimX=0; SimY=0; SimZ=0;    
    if dipoleaxis=='x'
        xpos_SPCW=xpos_SPCW+dx/2;
        xpos_dip=xpos_dip+dx/2;
        SimX=SimX+dx/2;        
    elseif dipoleaxis=='y'
        ypos_SPCW=ypos_SPCW+dy/2;
        ypos_dip=ypos_dip+dy/2;
        SimY=SimY+dy/2;        
    elseif dipoleaxis=='z'
        zpos_SPCW=zpos_SPCW+dz/2;
        zpos_dip=zpos_dip+dz/2;
        SimZ=SimZ+dz/2;
    end
    STLfile=[CurrentFolder '/2nmDisorder_SPCW_v5_STL/2nmDisorder_SPCW_v5_80unitcell.stl']
    lSTLimport(h,true,STLfile,nSiN,xpos_SPCW,ypos_SPCW,zpos_SPCW,toggle)
    %% Dipole source
    thetad=[90,90,0];
    phid=[90,0,0];
    DipoleSource(h,true,xpos_dip,ypos_dip,zpos_dip,freq,bandwidth,thetad(id),phid(id));    
    %% Setting FDTD simulation area
    AutoShutOff=1e-11;
    Courant=0.8;
    lFDTDcond(h,true,meshacc,timefactor,SimX,SimY,SimZ,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,dx,dy,dz,AutoShutOff,Courant);    
    %% Monitors
    xpos=0; ypos=0; zpos=0;
    lFDTDmonitors(h,true,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,xpos,ypos,zpos,freqhighres,timefactor)
    lindex(h,true,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,0,0,0)    
    %% Save and Run
    lcommand(h,['save("' SimName '_' dipoleaxis '")'],true);
    lcommand(h,'run',true);
    lcommand(h,['save("' SimName '_' dipoleaxis '")'],true);
    lcommand(h,'GFresult=getresult("GreenFunctionDATA","E")',true);
    lcommand(h,'n = getresult("index_monitor","index")',true);
    lcommand(h,'mu = getdata("DSource","moment")',true);
    lcommand(h,['matlabsave("muresult' dipoleaxis '",mu)'],true);
    lcommand(h,['matlabsave("nresult' num2str(dipoleaxis) '",n)'],true);
    lcommand(h,['matlabsave("GFresult' num2str(dipoleaxis) '",GFresult)'],true);
    %% Analysis
    targetx=0;targety=0;
    nresult=load(['nresult' num2str(dipoleaxis) '.mat']);
    nresult=nresult.n;
    muresult=load(['muresult' dipoleaxis '.mat']);
    muresult=muresult.mu;
    GFresult=load(['GFresult' num2str(dipoleaxis) '.mat']);
    GFresult=GFresult.GFresult;
    [Xplot,Yplot,GxN,GyN,GzN,zeropos]=Efieldplot(GFresult,muresult,fcs_d2,targetx,targety);
    [ffield,GxNf,GyNf,GzNf]=Efieldfrequencyplot(GFresult,zeropos,muresult);    
    lcommand(h,'switchtolayout',true);    
    toggle=2;
%     lgeometry(h,nSiN,dipoleaxis,dx,dy,dz,toggle)
    lSTLimport(h,true,STLfile,nSiN,xpos_SPCW,ypos_SPCW,zpos_SPCW,toggle)
    lcommand(h,['save("Vacuum' dipoleaxis '")'],true);
    lcommand(h,'run',true);
    lcommand(h,['save("Vacuum' dipoleaxis '")'],true);
    lcommand(h,'GFresult0=getresult("GreenFunctionDATA","E")',true);
    lcommand(h,['matlabsave("GFresult0' num2str(dipoleaxis) '",GFresult0)'],true);
    GFresult0=load(['GFresult0' num2str(dipoleaxis) '.mat']);
    GFresult0=GFresult0.GFresult0;
    [Xplot,Yplot,Gx0,Gy0,Gz0,zeropos]=Efieldplot(GFresult0,muresult,fcs_d2,targetx,targety);
    [ffield,Gx0f,Gy0f,Gz0f]=Efieldfrequencyplot(GFresult0,zeropos,muresult);
        
    if dipoleaxis=='x'
        GN(1,1,:)=GxN; GN(1,2,:)=GyN; GN(1,3,:)=conj(GzN);
        GNf(1,1,:)=GxNf; GNf(1,2,:)=GyNf; GNf(1,3,:)=conj(GzNf);
        G0(1,1,:)=Gx0; G0(1,2,:)=Gy0; G0(1,3,:)=Gz0;
        G0f(1,1,:)=Gx0f; G0f(1,2,:)=Gy0f; G0f(1,3,:)=Gz0f;
    elseif dipoleaxis=='y'
        GN(2,1,:)=conj(GxN); GN(2,2,:)=GyN; GN(2,3,:)=conj(GzN);
        GNf(2,1,:)=conj(GxNf); GNf(2,2,:)=GyNf; GNf(2,3,:)=conj(GzNf);
        G0(2,1,:)=Gx0; G0(2,2,:)=Gy0; G0(2,3,:)=Gz0;
        G0f(2,1,:)=Gx0f; G0f(2,2,:)=Gy0f; G0f(2,3,:)=Gz0f;
    elseif dipoleaxis=='z'
        GN(3,1,:)=GxN; GN(3,2,:)=GyN; GN(3,3,:)=GzN;
        GNf(3,1,:)=GxNf; GNf(3,2,:)=GyNf; GNf(3,3,:)=GzNf;
        G0(3,1,:)=Gx0; G0(3,2,:)=Gy0; G0(3,3,:)=Gz0;
        G0f(3,1,:)=Gx0f; G0f(3,2,:)=Gy0f; G0f(3,3,:)=Gz0f;
    end
    
    save('GreenN2D.mat','GN')
    save('GreenNf2D.mat','GNf')
    save('Green02D.mat','G0')
    save('Green0f2D.mat','G0f')    
    lcommand(h,'switchtolayout',true);
    appclose(h)
    
end
% cd(workspace)
% pltGF
% PltGammaJij
close all
