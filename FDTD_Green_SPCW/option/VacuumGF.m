function [G0,G0f]=VacuumGF(oneD,twoD,wSimVolx,wSimVoly,wSimVolz,dx,freq,bandwidth,freqhighres,timefactor)
close all
for id=1:3
    if id==1
        dipoleaxis='x'
    elseif id==2
        dipoleaxis='y'
    elseif id==3
        dipoleaxis='z'
    end
    fcs_d1=335.116e12;
    fcs_d2=351.726e12;
    cc=299792458; %speed of light
    %% Add Lumerical API to MATLAB
    CurrentFolder=pwd;
    path(path, 'C:\Program Files\Lumerical\FDTD\api\matlab')
    %appclose(h);
    h=appopen('fdtd');
    
    %% Load Lumerical Header API
    lheader(h,CurrentFolder)
    
    %% Simulation Parameters
    meshacc = 6;
    dy=dx;dz=dx;
    xSimMin=-wSimVolx;
    xSimMax=wSimVolx;
    ySimMin=-wSimVoly;
    ySimMax=wSimVoly;
    zSimMax=wSimVolz;
    zSimMin=-wSimVolz;
    
    %% #For dipole rotations
    thid=[90,90,0];
    phid=[0,90,0];
    
    %% Dipole Source
    thetad=[90,90,0];
    phid=[0,90,0];
    
    %% Create the geometry
    toggle=2;
    lgeometry(h,1,200e-9,dipoleaxis,dx,dy,dz,toggle)
    if dipoleaxis=='x'        
        DipoleSource(h,true,dx/2,0,0,freq,bandwidth,thetad(id),phid(id),dipoleaxis);
    elseif dipoleaxis=='y'
        DipoleSource(h,true,0,dy/2,0,freq,bandwidth,thetad(id),phid(id),dipoleaxis);
    elseif dipoleaxis=='z'
        DipoleSource(h,true,0,0,dz/2,freq,bandwidth,thetad(id),phid(id),dipoleaxis);
    end
    
    %% Setting FDTD simulation area
    Tf=timefactor*10e-6/cc;
    lFDTDcond(h,true,meshacc,Tf,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,dx,dy,dz);
%     lset(h,true,'x min bc','Symmetric')
%     lset(h,true,'y min bc','Anti-Symmetric')
%     lset(h,true,'z min bc','Symmetric')
    xpos=0;
    ypos=0;
    zpos=0;
    lFDTDmonitors(h,true,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,xpos,ypos,zpos,0,0,freqhighres,oneD,twoD)
    lindex(h,true,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,0,0,0)
    
    lcommand(h,['save("Vacuum' dipoleaxis '")'],true);
    lcommand(h,'run',true);
    lcommand(h,['save("Vacuum' dipoleaxis '")'],true);
    lcommand(h,'GFresult=getresult("GreenFunctionDATA","E")',true);
    lcommand(h,'delta_x=getdata("GreenFunctionDATA","delta_x")',true);
    lcommand(h,'delta_y=getdata("GreenFunctionDATA","delta_y")',true);
    lcommand(h,'n = getresult("index_monitor","index")',true);
    lcommand(h,'mu = getdata("DSource","moment")',true);
    lcommand(h,['matlabsave("muresult0' dipoleaxis '",mu)'],true);
    lcommand(h,['matlabsave("nresult0' num2str(dipoleaxis) '",n)'],true);
    lcommand(h,['matlabsave("GFresult0' num2str(dipoleaxis) '",GFresult)'],true);
    lcommand(h,['matlabsave("delta_x0' num2str(dipoleaxis) '",delta_x)'],true);
    lcommand(h,['matlabsave("delta_y0' num2str(dipoleaxis) '",delta_y)'],true);
    
    %% Analysis
    targetx=0;targety=0;
    nresult=load(['nresult0' num2str(dipoleaxis) '.mat']);
    nresult=nresult.n;
    muresult=load(['muresult0' dipoleaxis '.mat']);
    muresult=muresult.mu;
    delta_x=load(['delta_x0' dipoleaxis '.mat']);
    delta_x=delta_x.delta_x;
    delta_y=load(['delta_y0' dipoleaxis '.mat']);
    delta_y=delta_y.delta_y;
    GFresult=load(['GFresult0' num2str(dipoleaxis) '.mat']);
    GFresult=GFresult.GFresult;
    [Xplot,Yplot,Gx0,Gy0,Gz0,zeropos]=Efieldplot(GFresult,muresult,fcs_d2,targetx,targety,delta_x,delta_y,oneD,twoD,dy);
    [ffield,Gx0f,Gy0f,Gz0f]=Efieldfrequencyplot(GFresult,zeropos,muresult,oneD,twoD);
    
    if oneD
        if id==1
            G0(1,1)=Gx0; G0(1,2)=Gy0; G0(1,3)=Gz0;
            G0f(1,1,:)=Gx0f; G0f(1,2,:)=Gy0f; G0f(1,3,:)=Gz0f;
        elseif id==2
            G0(2,1)=Gx0; G0(2,2)=Gy0; G0(2,3)=Gz0;
            G0f(2,1,:)=Gx0f; G0f(2,2,:)=Gy0f; G0f(2,3,:)=Gz0f;
        elseif id==3
            G0(3,1)=Gx0; G0(3,2)=Gy0; G0(3,3)=Gz0;
            G0f(3,1,:)=Gx0f; G0f(3,2,:)=Gy0f; G0f(3,3,:)=Gz0f;
        end
        save('Green01D.mat','G0')
        save('Green0f1D.mat','G0f')
        
    elseif twoD
        if id==1
            G0(1,1,:)=Gx0; G0(1,2,:)=Gy0; G0(1,3,:)=Gz0;
            G0f(1,1,:)=Gx0f; G0f(1,2,:)=Gy0f; G0f(1,3,:)=Gz0f;
        elseif id==2
            G0(2,1,:)=Gx0; G0(2,2,:)=Gy0; G0(2,3,:)=Gz0;
            G0f(2,1,:)=Gx0f; G0f(2,2,:)=Gy0f; G0f(2,3,:)=Gz0f;
        elseif id==3
            G0(3,1,:)=Gx0; G0(3,2,:)=Gy0; G0(3,3,:)=Gz0;
            G0f(3,1,:)=Gx0f; G0f(3,2,:)=Gy0f; G0f(3,3,:)=Gz0f;
        end
        save('Green02D.mat','G0')
        save('Green0f2D.mat','G0f')
    end
    lcommand(h,'switchtolayout',true);
    appclose(h)
end