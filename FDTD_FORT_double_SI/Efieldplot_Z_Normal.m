function [Xplot,Yplot,Ex,Ey,Ez]=Efieldplot_Z_Normal(SILresult,nresult,ftrap)
close all
cc=299792458;eps0=8.854e-12;
Efield=SILresult.E;
ffield=SILresult.f;
nresult_x=nresult.index_x;
tmp=abs(ffield-ftrap);
[idx idx]=min(tmp);
closest=ffield(idx);
w=2*pi*closest;

[Xplot,Yplot]=meshgrid(SILresult.x,SILresult.y);
[Xsize,Ysize]=size(Xplot);
Ex=zeros(size(Xplot));
irun=0;

for i1=1:Xsize
    for i2=1:Ysize
        irun=irun+1;
        Ex(i1,i2)=squeeze((Efield(irun,1,idx)));
    end
end
ExNorm=(Ex-min(min(Ex)))/(max(max(Ex))-min(min(Ex)));

Ey=zeros(size(Xplot));
irun=0;
for i1=1:Xsize
    for i2=1:Ysize
        irun=irun+1;
        Ey(i1,i2)=squeeze((Efield(irun,2,idx)));
    end
end
EyNorm=(Ey-min(min(Ey)))/(max(max(Ey))-min(min(Ey)));

Ez=zeros(size(Xplot));
irun=0;
for i1=1:Xsize
    for i2=1:Ysize
        irun=irun+1;
        Ez(i1,i2)=squeeze((Efield(irun,3,idx)));
    end
end
EzNorm=(Ez-min(min(Ez)))/(max(max(Ez))-min(min(Ez)));

Nx=zeros(size(Xplot));
irun=0;
for i1=1:Xsize
    for i2=1:Ysize
        irun=irun+1;
        Nx(i1,i2)=squeeze((nresult_x(irun,1)));
    end
end

lambda=cc/ftrap*10^9;
A=figure(1)
ax1=axes;
cc=[0:0.01:1]
[C H1]=contourf(ax1,Xplot*10^6,Yplot*10^6,(abs(ExNorm).^2),cc);
set(H1,'LineColor','none')
colormap(ax1,jet)
caxis(ax1,[0 1])
hold on
[C H2]=contourf(Xplot*10^6,Yplot*10^6,Nx,[1.5 2])
hold off
hXLabel=xlabel(ax1,'y ($\mu$m)','Fontsize',15,'Interpreter','Latex')
hYLabel=ylabel(ax1,'z ($\mu$m)','Fontsize',15,'Interpreter','Latex')
set(H2,'color',[0.8 0.8 0.8],'Edgecolor',[0 0 0])
daspect(ax1,[1 1 1])
colorbar
ax1.XGrid='on'
ax1.YGrid='on'
ax1.FontSize=12;
title(['SIL APCW (' num2str(lambda,4) ' nm) X-pol'])
save(['Ex_XZ' num2str(lambda,4) 'nm.mat'],'Ex')
saveas(A,['Ix_XY' num2str(lambda,4) 'nm.png'])

B=figure(2)
ax1=axes;
cc=[0:0.01:1]
[C H1]=contourf(ax1,Xplot*10^6,Yplot*10^6,(abs(EyNorm).^2),cc);
set(H1,'LineColor','none')
colormap(ax1,jet)
caxis(ax1,[0 1])
hold on
[C H2]=contourf(Xplot*10^6,Yplot*10^6,Nx,[1.5 2])
hold off
hXLabel=xlabel(ax1,'y ($\mu$m)','Fontsize',15,'Interpreter','Latex')
hYLabel=ylabel(ax1,'z ($\mu$m)','Fontsize',15,'Interpreter','Latex')
set(H2,'color',[0.8 0.8 0.8],'Edgecolor',[0 0 0])
daspect(ax1,[1 1 1])
colorbar
ax1.XGrid='on'
ax1.YGrid='on'
ax1.FontSize=12;
title(['SIL APCW (' num2str(lambda,4) ' nm) Y-pol'])
save(['Ey_XZ' num2str(lambda,4) 'nm.mat'],'Ey')
saveas(B,['Iy_XY' num2str(lambda,4) 'nm.png'])

G=figure(3)
ax1=axes;
cc=[0:0.01:1]
[C H1]=contourf(ax1,Xplot*10^6,Yplot*10^6,(abs(EzNorm).^2),cc);
set(H1,'LineColor','none')
colormap(ax1,jet)
caxis(ax1,[0 1])
hold on
[C H2]=contourf(Xplot*10^6,Yplot*10^6,Nx,[1.5 2])
hold off
hXLabel=xlabel(ax1,'y ($\mu$m)','Fontsize',15,'Interpreter','Latex')
hYLabel=ylabel(ax1,'z ($\mu$m)','Fontsize',15,'Interpreter','Latex')
set(H2,'color',[0.8 0.8 0.8],'Edgecolor',[0 0 0])
daspect(ax1,[1 1 1])
colorbar
ax1.XGrid='on'
ax1.YGrid='on'
ax1.FontSize=12;
title(['SIL APCW (' num2str(lambda,4) ' nm) Z-pol'])
save(['Ez_XZ' num2str(lambda,4) 'nm.mat'],'Ez')
saveas(G,['Iz_XY' num2str(lambda,4) 'nm.png'])

end
