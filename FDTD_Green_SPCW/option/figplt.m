function figplt(Xplot,Yplot,ffield,prx,pry,GN,G0,GNf,G0f,GammaJijN,GammaJij0,oneD,twoD)
ex=[1 0 0];
ey=[0 1 0];
ez=[0 0 1];
etheta=(1/sqrt(2))*(ex+ez);
eeps=(1/sqrt(2))*(ex+1i*ez);
for i=1:size(GN,3)
    for j=1:size(GN,4)
        GammaJijNtheta(i,j)=etheta*squeeze(GammaJijN(:,:,i,j))*etheta';
        GammaJij0theta(i,j)=etheta*squeeze(GammaJij0(:,:,i,j))*etheta';
        GammaJijNeps(i,j)=eeps*squeeze(GammaJijN(:,:,i,j))*eeps';
        GammaJij0eps(i,j)=eeps*squeeze(GammaJij0(:,:,i,j))*eeps';
    end
end
if oneD
    figure(1)
    h1=plot(prx,squeeze(imag(GammaJijN(1,1,:,1))./imag(GammaJij0(1,1,:,1))));
    hTitle  = title ('Decay rate_{xx}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$Gamma1D/Gamma0_{xx}$$'],'Interpreter','latex');
    set(h1,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h1,['Gamma_xx.png'])
    
    figure(2)
    h2=plot(prx,squeeze(imag(GammaJijN(2,2,:,1))./imag(GammaJij0(2,2,:,1))));
    hTitle  = title ('Decay rate_{yy}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$Gamma1D/Gamma0_{yy}$$'],'Interpreter','latex');
    set(h2,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h2,['Gamma_yy.png'])
    
    figure(3)
    h3=plot(prx,squeeze(imag(GammaJijN(3,3,:,1))./imag(GammaJij0(3,3,:,1))));
    hTitle  = title ('Decay rate_{zz}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$Gamma1D/Gamma0_{zz}$$'],'Interpreter','latex');
    set(h3,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h3,['Gamma_zz.png'])
    
    figure(4)
    h4=plot(prx,squeeze(real(GammaJijN(1,1,:,1))-real(GammaJij0(1,1,:,1)))./(2*pi));
    hTitle  = title ('Lambshift_{xx}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$J1D-J0_{xx}$$'],'Interpreter','latex');
    set(h4,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h4,['Jij_xx.png'])
    
    figure(5)
    h5=plot(prx,squeeze(real(GammaJijN(2,2,:,1))-real(GammaJij0(2,2,:,1)))./(2*pi));
    hTitle  = title ('Lambshift_{yy}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$J1D-J0_{yy}$$'],'Interpreter','latex');
    set(h5,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h5,['Jij_yy.png'])
    
    figure(6)
    h6=plot(prx,squeeze(real(GammaJijN(3,3,:,1))-real(GammaJij0(3,3,:,1)))./(2*pi));
    hTitle  = title ('Lambshift_{zz}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$J1D-J0_{zz}$$'],'Interpreter','latex');
    set(h6,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h6,['Jij_zz.png'])
    
    figure(7)
    h7=plot(prx,squeeze(imag(GammaJijNeps(:,1))./imag(GammaJij0eps(:,1))));
    hTitle  = title ('Decay rate_{eps}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$Gamma_{total}/Gamma_{vac}$$'],'Interpreter','latex');
    set(h7,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h7,['Gamma_eps.png'])
    
    figure(8)
    h8=plot(prx,squeeze(imag(GammaJijNtheta(:,1))./imag(GammaJij0theta(:,1))));
    hTitle  = title ('Decay rate_{theta}');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$Gamma_{total}/Gamma_{vac}$$'],'Interpreter','latex');
    set(h8,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h8,['Gamma_theta.png'])
    
    figure(9)
    h9=plot(pry,squeeze(imag(GammaJijN(1,1,1,:))./imag(GammaJij0(1,1,1,:))));
    hTitle  = title ('Decay rate_{xx}');
    hXLabel = xlabel('y (m)');
    hYLabel = ylabel(['$$Gamma1D/Gamma0_{xx}$$'],'Interpreter','latex');
    set(h9,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h9,['Gamma_xx_yr.png'])
    
    figure(10)
    h10=plot(pry,squeeze(imag(GammaJijN(2,2,1,:))./imag(GammaJij0(2,2,1,:))));
    hTitle  = title ('Decay rate_{yy}');
    hXLabel = xlabel('y (m)');
    hYLabel = ylabel(['$$Gamma1D/Gamma0_{yy}$$'],'Interpreter','latex');
    set(h10,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h10,['Gamma_yy_yr.png'])
    
    figure(11)
    h11=plot(pry,squeeze(imag(GammaJijN(3,3,1,:))./imag(GammaJij0(3,3,1,:))));
    hTitle  = title ('Decay rate_{zz}');
    hXLabel = xlabel('y (m)');
    hYLabel = ylabel(['$$Gamma1D/Gamma0_{zz}$$'],'Interpreter','latex');
    set(h11,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h11,['Gamma_zz_yr.png'])

    figure(12)
    h12=plot(pry,squeeze(real(GammaJijN(1,1,1,:))-real(GammaJij0(1,1,1,:)))./(2*pi));
    hTitle  = title ('Lambshift_{xx}');
    hXLabel = xlabel('y (m)');
    hYLabel = ylabel(['$$J1D-J0_{xx}$$'],'Interpreter','latex');
    set(h12,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h12,['Jij_xx_yr.png'])
    
    figure(13)
    h13=plot(pry,squeeze(real(GammaJijN(2,2,1,:))-real(GammaJij0(2,2,1,:)))./(2*pi));
    hTitle  = title ('Lambshift_{yy}');
    hXLabel = xlabel('y (m)');
    hYLabel = ylabel(['$$J1D-J0_{yy}$$'],'Interpreter','latex');
    set(h13,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h13,['Jij_yy_yr.png'])
    
    figure(14)
    h14=plot(pry,squeeze(real(GammaJijN(3,3,1,:))-real(GammaJij0(3,3,1,:)))./(2*pi));
    hTitle  = title ('Lambshift_{zz}');
    hXLabel = xlabel('y (m)');
    hYLabel = ylabel(['$$J1D-J0_{zz}$$'],'Interpreter','latex');
    set(h14,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h14,['Jij_zz_yr.png'])
    
    
    px=[-1*fliplr(prx),prx];
    py=[-1*fliplr(pry),pry];
    Decay_xx=[rot90(squeeze(imag(GammaJijN(1,1,:,:))./imag(GammaJij0(1,1,:,:))),2), ...
        flipud(squeeze(imag(GammaJijN(1,1,:,:))./imag(GammaJij0(1,1,:,:)))); ...
        fliplr(squeeze(imag(GammaJijN(1,1,:,:))./imag(GammaJij0(1,1,:,:)))),squeeze(imag(GammaJijN(1,1,:,:))./imag(GammaJij0(1,1,:,:)))];
    Decay_yy=[rot90(squeeze(imag(GammaJijN(2,2,:,:))./imag(GammaJij0(2,2,:,:))),2), ...
        flipud(squeeze(imag(GammaJijN(2,2,:,:))./imag(GammaJij0(2,2,:,:)))); ...
        fliplr(squeeze(imag(GammaJijN(2,2,:,:))./imag(GammaJij0(2,2,:,:)))),squeeze(imag(GammaJijN(2,2,:,:))./imag(GammaJij0(2,2,:,:)))];
    Decay_zz=[rot90(squeeze(imag(GammaJijN(3,3,:,:))./imag(GammaJij0(3,3,:,:))),2), ...
        flipud(squeeze(imag(GammaJijN(3,3,:,:))./imag(GammaJij0(3,3,:,:)))); ...
        fliplr(squeeze(imag(GammaJijN(3,3,:,:))./imag(GammaJij0(3,3,:,:)))),squeeze(imag(GammaJijN(3,3,:,:))./imag(GammaJij0(3,3,:,:)))];
    Lambshift_xx=[rot90(squeeze(real(GammaJijN(1,1,:,:))-real(GammaJij0(1,1,:,:)))./(2*pi),2), ...
        flipud(squeeze(real(GammaJijN(1,1,:,:))-real(GammaJij0(1,1,:,:)))./(2*pi)); ...
        fliplr(squeeze(real(GammaJijN(1,1,:,:))-real(GammaJij0(1,1,:,:)))./(2*pi)),squeeze(real(GammaJijN(1,1,:,:))-real(GammaJij0(1,1,:,:)))./(2*pi)];
    Lambshift_yy=[rot90(squeeze(real(GammaJijN(2,2,:,:))-real(GammaJij0(2,2,:,:)))./(2*pi),2), ...
        flipud(squeeze(real(GammaJijN(2,2,:,:))-real(GammaJij0(2,2,:,:)))./(2*pi)); ...
        fliplr(squeeze(real(GammaJijN(2,2,:,:))-real(GammaJij0(2,2,:,:)))./(2*pi)),squeeze(real(GammaJijN(2,2,:,:))-real(GammaJij0(2,2,:,:)))./(2*pi)];
    Lambshift_zz=[rot90(squeeze(real(GammaJijN(3,3,:,:))-real(GammaJij0(3,3,:,:)))./(2*pi),2), ...
        flipud(squeeze(real(GammaJijN(3,3,:,:))-real(GammaJij0(3,3,:,:)))./(2*pi)); ...
        fliplr(squeeze(real(GammaJijN(3,3,:,:))-real(GammaJij0(3,3,:,:)))./(2*pi)),squeeze(real(GammaJijN(3,3,:,:))-real(GammaJij0(3,3,:,:)))./(2*pi)];
    
    figure(15)
    h15=surf(py,px,Decay_xx)
    title('Decay rate_{xx}')
    colorbar
    xlabel('y (m)')
    ylabel('x (m)')
    zlabel('Gamma_{total}/Gamma_{rad}')
    saveas(h15,'Decay_rate2D_xx.png')
    
    figure(16)
    h16=surf(py,px,Decay_yy)
    title('Decay rate_{yy}')
    colorbar
    xlabel('y (m)')
    ylabel('x (m)')
    zlabel('Gamma_{total}/Gamma_{rad}')
    saveas(h16,'Decay_rate2D_yy.png')

    figure(17)
    h17=surf(py,px,Decay_zz)
    title('Decay rate_{zz}')
    colorbar
    xlabel('y (m)')
    ylabel('x (m)')
    zlabel('Gamma_{total}/Gamma_{rad}')
    saveas(h17,'Decay_rate2D_zz.png')
    
    figure(18)
    h18=surf(py,px,Lambshift_xx)
    title('Lambshift_{xx}')
    colorbar
    xlabel('y (m)')
    ylabel('x (m)')
    zlabel('Lambshift_{xx} (Hz)')
    saveas(h18,'Lambshift2D_xx.png')
    
    figure(19)
    h19=surf(py,px,Lambshift_yy)
    title('Lambshift_{yy}')
    colorbar
    xlabel('y (m)')
    ylabel('x (m)')
    zlabel('Lambshift_{yy} (Hz)')
    saveas(h19,'Lambshift2D_yy.png')
    
    figure(20)
    h20=surf(py,px,Lambshift_zz)
    title('Lambshift_{zz}')
    colorbar
    xlabel('y (m)')
    ylabel('x (m)')
    zlabel('Lambshift_{zz} (Hz)')
    saveas(h20,'Lambshift2D_zz.png')
    
elseif twoD
    figure(1)
    %    [C,h1] =contourf(Xplot,Yplot,squeeze(imag(GN(1,1,:,:))),20);
    h1=plot(Xplot,squeeze(imag(GN(1,1,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{xx})$$'],'Interpreter','latex');
    set(h1,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h1,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h1,['Img_GF_xx_rrp.png'])
    
    figure(2)
    %    [C,h2] =contourf(Xplot,Yplot,squeeze(imag(GN(1,2,:,:))),20);
    h2=plot(Xplot,squeeze(imag(GN(1,2,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{xy})$$'],'Interpreter','latex');
    set(h2,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h2,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h2,['Img_GF_xy_rrp.png'])
    
    figure(3)
    %    [C,h3] = contourf(Xplot,Yplot,squeeze(imag(GN(1,3,:,:))),20);
    h3=plot(Xplot,squeeze(imag(GN(1,3,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{xz})$$'],'Interpreter','latex');
    set(h3,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h3,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h3,['Img_GF_xz_rrp.png'])
    
    figure(4)
    %    [C,h4] =contourf(Xplot,Yplot,squeeze(imag(GN(2,1,:,:))),20);
    h4=plot(Xplot,squeeze(imag(GN(2,1,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{yx})$$'],'Interpreter','latex');
    set(h4,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h4,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h4,['Img_GF_yx_rrp.png'])
    
    figure(5)
    %    [C,h5] =contourf(Xplot,Yplot,squeeze(imag(GN(2,2,:,:))),20);
    h5=plot(Xplot,squeeze(imag(GN(2,2,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{yy})$$'],'Interpreter','latex');
    set(h5,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h5,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h5,['Img_GF_yy_rrp.png'])
    
    figure(6)
    %    [C,h6] = contourf(Xplot,Yplot,squeeze(imag(GN(2,3,:,:))),20);
    h6=plot(Xplot,squeeze(imag(GN(2,3,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{yz})$$'],'Interpreter','latex');
    set(h6,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h6,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h6,['Img_GF_yz_rrp.png'])
    
    figure(7)
    %    [C,h7] = contourf(Xplot,Yplot,squeeze(imag(GN(3,1,:,:))),20);
    h7=plot(Xplot,squeeze(imag(GN(3,1,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{zx})$$'],'Interpreter','latex');
    set(h7,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h7,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h7,['Img_GF_zx_rrp.png'])
    
    figure(8)
    %    [C,h8] = contourf(Xplot,Yplot,squeeze(imag(GN(3,2,:,:))),20);
    h8=plot(Xplot,squeeze(imag(GN(2,1,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{zy})$$'],'Interpreter','latex');
    set(h8,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h8,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h8,['Img_GF_zy_rrp.png'])
    
    figure(9)
    %    [C,h9] = contourf(Xplot,Yplot,squeeze(imag(GN(3,3,:,:))),20);
    h9=plot(Xplot,squeeze(imag(GN(3,3,:,:))));
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$imag(Green_{zz})$$'],'Interpreter','latex');
    set(h9,'Color',[0 0 0],'LineWidth',1.5 );
    %    colorbar
    grid on
    %    daspect([1 1 1])
    %    set(h9,'LineColor','none')
    set(gca,'FontSize',5)
    saveas(h9,['Img_GF_zz_rrp.png'])
    
    figure(10)
    h10=plot(ffield,imag(squeeze(GNf(1,1,:)))./imag(squeeze(G0f(1,1,:))));
    hTitle  = title ('Normalized Green function');
    hXLabel = xlabel('frequency (Hz)');
    hYLabel = ylabel(['$$Im(G_{xx})/Im(G_{xx})_{vac}$$'],'Interpreter','latex');
    set(h10,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h10,['Img_GF_xx.png'])
    
    figure(11)
    h11=plot(ffield,imag(squeeze(GNf(2,2,:)))./imag(squeeze(G0f(2,2,:))));
    hTitle  = title ('Normalized Green function');
    hXLabel = xlabel('frequency (Hz)');
    hYLabel = ylabel(['$$Im(G_{yy})/Im(G_{yy})_{vac}$$'],'Interpreter','latex');
    set(h11,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h11,['Img_GF_yy.png'])
    
    figure(12)
    h12=plot(ffield,imag(squeeze(GNf(3,3,:)))./imag(squeeze(G0f(3,3,:))));
    hTitle  = title ('Normalized Green function');
    hXLabel = xlabel('frequency (Hz)');
    hYLabel = ylabel(['$$Im(G_{zz})/Im(G_{zz})_{vac}$$'],'Interpreter','latex');
    set(h12,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h12,['Img_GF_zz.png'])
    
    figure(13)
    h13=plot(ffield,imag(squeeze(GNf(1,2,:)))./imag(squeeze(G0f(1,2,:))));
    hTitle  = title ('Normalized Green function');
    hXLabel = xlabel('frequency (Hz)');
    hYLabel = ylabel(['$$Im(G_{xy})/Im(G_{xy})_{vac}$$'],'Interpreter','latex');
    set(h13,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h13,['Img_GF_xy.png'])
    
    figure(14)
    h14=plot(ffield,imag(squeeze(GNf(1,3,:)))./imag(squeeze(G0f(1,3,:))));
    hTitle  = title ('Normalized Green function');
    hXLabel = xlabel('frequency (Hz)');
    hYLabel = ylabel(['$$Im(G_{xz})/Im(G_{xz})_{vac}$$'],'Interpreter','latex');
    set(h14,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h14,['Img_GF_xz.png'])
    
    figure(15)
    h15=plot(ffield,imag(squeeze(GNf(2,3,:)))./imag(squeeze(G0f(2,3,:))));
    hTitle  = title ('Normalized Green function');
    hXLabel = xlabel('frequency (Hz)');
    hYLabel = ylabel(['$$Im(G_{yz})/Im(G_{yz})_{vac}$$'],'Interpreter','latex');
    set(h15,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h15,['Img_GF_xx.png'])
    
    figure(16)
    h16=plot(Xplot,squeeze(real(GN(1,1,:,:))-real(G0(1,1,:,:))));
    hTitle  = title ('real(Green function)');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$real(Green_{xx})$$'],'Interpreter','latex');
    set(h16,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h16,['real_GF_xx.png'])
    
    figure(17)
    h17=plot(Xplot,squeeze(real(GN(1,2,:,:))-real(G0(1,2,:,:))));
    hTitle  = title ('real(Green function)');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$real(Green_{xy})$$'],'Interpreter','latex');
    set(h17,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h17,['real_GF_xy.png'])
    
    figure(18)
    h18=plot(Xplot,squeeze(real(GN(1,3,:,:))-real(G0(1,3,:,:))));
    hTitle  = title ('real(Green function)');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$real(Green_{xz})$$'],'Interpreter','latex');
    set(h18,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h18,['real_GF_xz.png'])
    
    figure(19)
    h19=plot(Xplot,squeeze(real(GN(2,2,:,:))-real(G0(2,2,:,:))));
    hTitle  = title ('real(Green function)');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$real(Green_{yy})$$'],'Interpreter','latex');
    set(h19,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h19,['real_GF_yy.png'])
    
    figure(20)
    h20=plot(Xplot,squeeze(real(GN(2,3,:,:))-real(G0(2,3,:,:))));
    hTitle  = title ('real(Green function)');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$real(Green_{yz})$$'],'Interpreter','latex');
    set(h20,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h20,['real_GF_yz.png'])
    
    figure(21)
    h21=plot(Xplot,squeeze(real(GN(3,3,:,:))-real(G0(3,3,:,:))));
    hTitle  = title ('real(Green function)');
    hXLabel = xlabel('x (m)');
    hYLabel = ylabel(['$$real(Green_{zz})$$'],'Interpreter','latex');
    set(h21,'Color',[0 0 0],'LineWidth',1.5 );
    saveas(h21,['real_GF_zz.png'])
    
end


