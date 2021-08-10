function [Aeff_even,Aeff_odd]=Modearea(fname,a0,nband)
%% Effective mode area calculation (light propagation direction kx) %%
%% Constant
SiN=4;
c=299792458;
ep0=8.85e-12;
hbar=1.05457148e-34;
h=6.626e-34;
%  a0=374e-9;
%  nband=20;
%% Data load and Unit conversion
load('data.mat');
data=data{1};
SIunit_omega=(c/a0)*2*pi;  % Angular frequency in Hz
SIunit_energy=SIunit_omega*hbar;   % Energy in J
SIunit_prop=(2*pi/a0);
SIunit_mass=hbar/(a0^2*2*pi*c/a0);
SIunit_k0=0.5*2*pi/a0;
kpoints=data.vectors(:,1);
lightcone=data.vectors(:,1);

%% Plot preparation
% Xplot=linspace(data.size(1)*(-a0/2)*10^9,data.size(1)*(a0/2)*10^9,data.grid(1)); % unit (nm)
Yplot=linspace((data.size(2)*(-a0/2)*10^6)+(((-a0/2)*10^6)/data.grid(1)),(data.size(2)*(a0/2)*10^6)+(((-a0/2)*10^6)/data.grid(1)),data.grid(2)); % unit (nm)
Zplot=linspace((data.size(3)*(-a0/2)*10^6)+(((-a0/2)*10^6)/data.grid(1)),(data.size(3)*(a0/2)*10^6)+(((-a0/2)*10^6)/data.grid(1)),data.grid(3)); % unit (nm)
Xplot=linspace(0,data.size(1)*(a0)*10^6,data.grid(1)); % unit (nm)
%Yplot=linspace(0,data.size(2)*(a0)*10^9,data.grid(2)); % unit (nm)
%Zplot=linspace(0,data.size(3)*(a0)*10^9,data.grid(3)); % unit (nm)
dx=data.size(1)*a0/(data.grid(1));
dy=data.size(2)*a0/(data.grid(2));
dz=data.size(3)*a0/(data.grid(3));

%% Load electric field data over the volume of interest
kx_edge=length(kpoints);
Epsilon=h5read(strcat(fname,'-epsilon.h5'),'//data');


if nband<10
    if kx_edge<10
        Efield_x_r=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b0',num2str(nband),'.x.zeven.h5'),'//x.r');
        Efield_x_i=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b0',num2str(nband),'.x.zeven.h5'),'//x.i');
        Efield_y_r=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b0',num2str(nband),'.y.zeven.h5'),'//y.r');
        Efield_y_i=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b0',num2str(nband),'.y.zeven.h5'),'//y.i');
        Efield_z_r=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b0',num2str(nband),'.z.zeven.h5'),'//z.r');
        Efield_z_i=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b0',num2str(nband),'.z.zeven.h5'),'//z.i');
    else
        Efield_x_r=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b0',num2str(nband),'.x.zeven.h5'),'//x.r');
        Efield_x_i=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b0',num2str(nband),'.x.zeven.h5'),'//x.i');
        Efield_y_r=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b0',num2str(nband),'.y.zeven.h5'),'//y.r');
        Efield_y_i=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b0',num2str(nband),'.y.zeven.h5'),'//y.i');
        Efield_z_r=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b0',num2str(nband),'.z.zeven.h5'),'//z.r');
        Efield_z_i=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b0',num2str(nband),'.z.zeven.h5'),'//z.i');
    end
else nband>10
    if kx_edge<10
        Efield_x_r=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b',num2str(nband),'.x.zeven.h5'),'//x.r');
        Efield_x_i=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b',num2str(nband),'.x.zeven.h5'),'//x.i');
        Efield_y_r=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b',num2str(nband),'.y.zeven.h5'),'//y.r');
        Efield_y_i=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b',num2str(nband),'.y.zeven.h5'),'//y.i');
        Efield_z_r=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b',num2str(nband),'.z.zeven.h5'),'//z.r');
        Efield_z_i=h5read(strcat(fname,'-e.k0',num2str(kx_edge),'.b',num2str(nband),'.z.zeven.h5'),'//z.i');
    else
        Efield_x_r=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b',num2str(nband),'.x.zeven.h5'),'//x.r');
        Efield_x_i=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b',num2str(nband),'.x.zeven.h5'),'//x.i');
        Efield_y_r=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b',num2str(nband),'.y.zeven.h5'),'//y.r');
        Efield_y_i=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b',num2str(nband),'.y.zeven.h5'),'//y.i');
        Efield_z_r=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b',num2str(nband),'.z.zeven.h5'),'//z.r');
        Efield_z_i=h5read(strcat(fname,'-e.k',num2str(kx_edge),'.b',num2str(nband),'.z.zeven.h5'),'//z.i');
    end
end
Efield_x=Efield_x_r+1i*Efield_x_i;
Efield_y=Efield_y_r+1i*Efield_y_i;
Efield_z=Efield_z_r+1i*Efield_z_i;
clear Efield_x_i Efield_x_r Efield_y_i Efield_y_r Efield_z_i Efield_z_r

%% Effective mode area calculation
Int=abs(Efield_x).^2+abs(Efield_y).^2+abs(Efield_z).^2;
ModeAeff=sum(sum(sum(Epsilon.*Int.*dx.*dy.*dz)))./(a0.*Epsilon.*Int).*10^12;

%% Linear Interporlation
[My YmaxIndex]=min(abs(Yplot-0.22));
[Mz ZmaxIndex]=min(abs(Zplot-0.22));
[My YminIndex]=min(abs(Yplot+0.22));
[Mz ZminIndex]=min(abs(Zplot+0.22));
Yplt=Yplot(YminIndex:YmaxIndex);
Zplt=Zplot(ZminIndex:ZmaxIndex);
Xplt=Xplot;
AreaYZplotE=squeeze(ModeAeff(ZminIndex:ZmaxIndex,length(Xplot)/2,YminIndex:YmaxIndex));
AreaYZplotO=squeeze(ModeAeff(ZminIndex:ZmaxIndex,1,YminIndex:YmaxIndex));
AreaXYplot=squeeze(ModeAeff(length(Zplot)/2,:,YminIndex:YmaxIndex));
EpsYZplot=squeeze(Epsilon(ZminIndex:ZmaxIndex,length(Xplot)/2,YminIndex:YmaxIndex));
EpsXYplot=squeeze(Epsilon(length(Zplot)/2,:,YminIndex:YmaxIndex));

[ZZ,YY]=meshgrid(Zplt,Yplt);
[ZZ1,YY1]=meshgrid(Zplt(1):0.01:Zplt(end),Yplt(1):0.01:Yplt(end));
[XXt,YYt]=meshgrid(Xplt,Yplt);
[XXt1,YYt1]=meshgrid(Xplt(1):0.01:Xplt(end),Yplt(1):0.01:Yplt(end));
YYplt=linspace(Yplt(1),Yplt(end),size(YY1,1));
ZZplt=linspace(Zplt(1),Zplt(end),size(ZZ1,2));
Aplt_even=interp2(ZZ,YY,AreaYZplotE',ZZ1,YY1,'linear');
Aplt_odd=interp2(ZZ,YY,AreaYZplotO',ZZ1,YY1,'linear');
AAplt=interp2(XXt,YYt,AreaXYplot',XXt1,YYt1,'linear');
EpsXY=interp2(XXt,YYt,EpsXYplot',XXt1,YYt1,'linear');
EpsYZ=interp2(ZZ,YY,EpsYZplot',ZZ1,YY1,'linear');


XXplt=linspace(Xplt(1),Xplt(end),size(XXt1,2));
YYplt=linspace(Yplt(1),Yplt(end),size(YYt1,1));
XXplot=repmat([-fliplr(XXplt),XXplt],length(YYplt),1);
YYplot=repmat(YYplt,length(XXplt)*2,1);
AArea=[fliplr(AAplt),AAplt];
EEps=[fliplr(EpsXY),EpsXY];

clear caxis
%% Plot
A=figure(4)
cc=[0:0.005:1.1];
zlevs=[0.18 0.2];
zlevs1=[0.3 0.6];
[Y,Z]=meshgrid(Yplot,Zplot);
ax3=axes;
[C,h25]=contourf(ax3,YY1',ZZ1',Aplt_even',cc);
set(h25,'LineColor','none')
caxis([0 1.1])
% shading interp
hold on
[C,h]=contour(ax3,YY1',ZZ1',Aplt_even',zlevs,'-','color',[0 0 0],'linewidth',1.5)
[C,h]=contour(ax3,YY1',ZZ1',Aplt_even',zlevs1,'-','color',[0 0 0],'linewidth',1.5)
text(-0.015,0.053,['0.2'],'FontSize',12,'Color',[0 0 0])
text(-0.08,0,['0.18'],'FontSize',12,'Color',[0 0 0])
text(-0.015,0.12,['0.3'],'FontSize',12)
text(-0.015,0.165,['0.6'],'FontSize',12)
[C,h26]=contourf(YY1',ZZ1',EpsYZ',[3 4.5])
hold off
hXLabel = xlabel(ax3,['y $(\mu m)$'],'Fontsize',15,'Interpreter','Latex');
hYLabel = ylabel(ax3,['z $(\mu m)$'],'Fontsize',15,'Interpreter','Latex');
xlim(ax3,[YY1(1),YY1(end)])
ylim(ax3,[ZZ1(1),ZZ1(end)])
colormap(ax3,flipud(hot))
set(h26,'color',[0.8 0.8 0.8],'Edgecolor',[0 0 0])
daspect(ax3,[1 1 1])
h1=colorbar
ylabel(h1,['$$A_{\rm eff} (\mu m^2)$$'],'Rotation',270,'Fontsize',14,'Interpreter','Latex')
hh=get(h1,'ylabel');
oldpos = get(hh, 'Position');
set(hh, 'Position', oldpos + [1.8, 0, 0])
ax3.XGrid='on'
ax3.YGrid='on'
ax3.XTick=[-0.2:0.1:0.2];
ax3.YTick=[-0.2:0.1:0.2];
ax3.FontSize=12
set(gca,'FontSize',11)
set(gca,'Fontname','times')
grid on
saveas(A,['Modearea_even' num2str(nband) '.png'])

B=figure(5)
cc=[0:0.005:1.1];
zlevs=[0.18 0.2];
zlevs1=[0.3 0.6];
ax3=axes;
[C,h25]=contourf(ax3,YY1',ZZ1',Aplt_odd',cc);
set(h25,'LineColor','none')
caxis([0 1.1])
hold on
[C,h]=contour(ax3,YY1',ZZ1',Aplt_odd',zlevs,'-','color',[0 0 0],'linewidth',1.5)
[C,h]=contour(ax3,YY1',ZZ1',Aplt_odd',zlevs1,'-','color',[0 0 0],'linewidth',1.5)
text(-0.015,0.053,['0.2'],'FontSize',12,'Color',[0 0 0])
text(-0.08,0,['0.18'],'FontSize',12,'Color',[0 0 0])
text(-0.015,0.12,['0.3'],'FontSize',12)
text(-0.015,0.165,['0.6'],'FontSize',12)
[C,h26]=contourf(YY1',ZZ1',EpsYZ',[3 4.5])
hold off
hXLabel = xlabel(ax3,['y $(\mu m)$'],'Fontsize',15,'Interpreter','Latex');
hYLabel = ylabel(ax3,['z $(\mu m)$'],'Fontsize',15,'Interpreter','Latex');
xlim(ax3,[YY1(1),YY1(end)])
ylim(ax3,[ZZ1(1),ZZ1(end)])
colormap(ax3,flipud(hot))
set(h26,'color',[0.8 0.8 0.8],'Edgecolor',[0 0 0])
daspect(ax3,[1 1 1])
h1=colorbar
ylabel(h1,['$$A_{\rm eff} (\mu m^2)$$'],'Rotation',270,'Fontsize',14,'Interpreter','Latex')
hh=get(h1,'ylabel');
oldpos = get(hh, 'Position');
set(hh, 'Position', oldpos + [1.8, 0, 0])
ax3.XGrid='on'
ax3.YGrid='on'
ax3.XTick=[-0.2:0.1:0.2];
ax3.YTick=[-0.2:0.1:0.2];
ax3.FontSize=12
set(gca,'FontSize',11)
set(gca,'Fontname','times')
grid on
saveas(B,['Modearea_odd' num2str(nband) '.png'])

C=figure(6)
cc=[0:0.005:1.1];
zlevs=[0.18 0.2];
zlevs1=[0.3 0.6];
ax3=axes;
[C,h25]=contourf(ax3,XXplot,YYplot',AArea,cc);
set(h25,'LineColor','none')
caxis([0 1.1])
hold on
[C,h]=contour(ax3,XXplot,YYplot',AArea,zlevs,'-','color',[0 0 0],'linewidth',1.5)
[C,h]=contour(ax3,XXplot,YYplot',AArea,zlevs1,'-','color',[0 0 0],'linewidth',1.5)
text(-0.015,0.053,['0.2'],'FontSize',12,'Color',[0 0 0])
text(-0.08,0,['0.18'],'FontSize',12,'Color',[0 0 0])
text(-0.015,0.12,['0.3'],'FontSize',12)
text(-0.015,0.165,['0.6'],'FontSize',12)
[C,h26]=contourf(XXplot,YYplot',EEps,[3 4.5])
hold off
hXLabel = xlabel(ax3,['x $(\mu m)$'],'Fontsize',15,'Interpreter','Latex');
hYLabel = ylabel(ax3,['y $(\mu m)$'],'Fontsize',15,'Interpreter','Latex');
xlim(ax3,[XXplot(1),XXplot(end)])
ylim(ax3,[YYplot(1),YYplot(end)])
colormap(ax3,flipud(hot))
set(h26,'color',[0.8 0.8 0.8],'Edgecolor',[0 0 0])
daspect(ax3,[1 1 1])
h1=colorbar
ylabel(h1,['$$A_{\rm eff} (\mu m^2)$$'],'Rotation',270,'Fontsize',14,'Interpreter','Latex')
hh=get(h1,'ylabel');
oldpos = get(hh, 'Position');
set(hh, 'Position', oldpos + [1.8, 0, 0])
ax3.XGrid='on'
ax3.YGrid='on'
ax3.XTick=[-0.2:0.1:0.2];
ax3.YTick=[-0.2:0.1:0.2];
ax3.FontSize=12
set(gca,'FontSize',11)
set(gca,'Fontname','times')
grid on

Aeff_even=ModeAeff(round(size(ModeAeff,1)/2),1,round(size(ModeAeff,3)/2))
Aeff_odd=ModeAeff(round(size(ModeAeff,1)/2),round(size(ModeAeff,2)/2),round(size(ModeAeff,3)/2))

