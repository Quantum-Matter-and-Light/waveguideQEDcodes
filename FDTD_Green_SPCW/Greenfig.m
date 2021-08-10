function figplt(Xplot,Yplot,ffield,GN,G0,GNf,G0f)

figure(1)
h1=plot(ffield,imag(squeeze(GNf(1,1,:)))./imag(squeeze(G0f(1,1,:))));
hTitle  = title ('Normalized Green function','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Im(G_{xx})/Im(G_{xx})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(h1,'Color',[0 0 0],'LineWidth',1.5 );
set(gca,'Fontsize',12)
saveas(h1,['Img_GF_xx.png'])

figure(2)
h2=plot(ffield,imag(squeeze(GNf(2,2,:)))./imag(squeeze(G0f(2,2,:))));
hTitle  = title ('Normalized Green function','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Im(G_{yy})/Im(G_{yy})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h2,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h2,['Img_GF_yy.png'])

figure(3)
h3=plot(ffield,imag(squeeze(GNf(3,3,:)))./imag(squeeze(G0f(3,3,:))));
hTitle  = title ('Normalized Green function','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Im(G_{zz})/Im(G_{zz})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h3,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h3,['Img_GF_zz.png'])

figure(4)
h4=plot(ffield,imag(squeeze(GNf(1,2,:)))./imag(squeeze(G0f(1,2,:))));
hTitle  = title ('Normalized Green function','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Im(G_{xy})/Im(G_{xy})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h4,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h4,['Img_GF_xy.png'])

figure(5)
h5=plot(ffield,imag(squeeze(GNf(1,3,:)))./imag(squeeze(G0f(1,3,:))));
hTitle  = title ('Normalized Green function','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Im(G_{xz})/Im(G_{xz})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h5,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h5,['Img_GF_xz.png'])

figure(6)
h6=plot(ffield,imag(squeeze(GNf(2,3,:)))./imag(squeeze(G0f(2,3,:))));
hTitle  = title ('Normalized Green function','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Im(G_{yz})/Im(G_{yz})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h6,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h6,['Img_GF_yz.png'])

figure(7)
h7=plot(ffield,squeeze(real(GNf(1,1,:,:))-real(G0f(1,1,:,:))));
hTitle  = title ('Re[Inhomogeneous Green function]','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Re(Green_{xx})-Re(Green_{xx})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h7,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h7,['real_GF_xx.png'])

figure(8)
h8=plot(ffield,squeeze(real(GNf(1,2,:,:))-real(G0f(1,2,:,:))));
hTitle  = title ('Re[Inhomogeneous Green function]','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Re(Green_{xy})-Re(Green_{xy})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h8,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h8,['real_GF_xy.png'])

figure(9)
h9=plot(ffield,squeeze(real(GNf(1,3,:,:))-real(G0f(1,3,:,:))));
hTitle  = title ('Re[Inhomogeneous Green function]','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Re(Green_{xz})-Re(Green_{xz})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h9,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h9,['real_GF_xz.png'])

figure(10)
h10=plot(ffield,squeeze(real(GNf(2,2,:,:))-real(G0f(2,2,:,:))));
hTitle  = title ('Re[Inhomogeneous Green function]','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Re(Green_{yy})-Re(Green_{yy})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h10,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h10,['real_GF_yy.png'])

figure(11)
h11=plot(ffield,squeeze(real(GNf(2,3,:,:))-real(G0f(2,3,:,:))));
hTitle  = title ('Re[Inhomogeneous Green function]','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Re(Green_{yz})-Re(Green_{yz})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h11,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h11,['real_GF_yz.png'])

figure(12)
h12=plot(ffield,squeeze(real(GNf(3,3,:,:))-real(G0f(3,3,:,:))));
hTitle  = title ('Re[Inhomogeneous Green function]','Fontsize',14);
hXLabel = xlabel('frequency (Hz)','Interpreter','latex','Fontsize',14);
hYLabel = ylabel(['$$Re(Green_{zz})-Re(Green_{zz})_{vac}$$'],'Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12)
set(h12,'Color',[0 0 0],'LineWidth',1.5 );
saveas(h12,['real_GF_zz.png'])

end