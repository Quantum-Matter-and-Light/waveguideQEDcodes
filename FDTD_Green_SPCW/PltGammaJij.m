clear all
close all
%% Load file
load('Gamma_tot.mat')
load('Gamma_vac.mat')
load('Jij.mat')
load('Jij_vac.mat')
load('Gamma1D.mat')
load('freq.mat')
load('Xplot.mat')
load('NLoJij0.mat')
load('NLoJij.mat')
load('NLoGamma.mat')
load('NLoGamma0.mat')
%% Plot
[M I]=max(squeeze(Gamma_yy));
fg=ffield(I);
for i=5:35
    [XX(i) Index(i-4)]=min(abs(Xplot+Xplot(end)-(370e-9)*1.0048*i));
end
%% Gamma and Jij wide range plot
A=figure(1)
set(A,'PaperPositionMode','auto')
set(A,'Position',[1000 800 1400 1000]);
H=subplot(4,2,2,'Position',[0.7 0.76 0.22 0.21])
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma0),'k-','linewidth',1.5)
hold on
hold off
xlim([-0.1 1])
ylim([1 5000])
grid on
set(gca,'Fontsize',12)
set(gca,'xticklabel',{[]}) 
set(gca,'YMinorTick','on')
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')

H=subplot(4,2,1,'Position',[0.08 0.76 0.58 0.21])
ha = area([384.1 ffield(end)*10^-12], [5000 5000]); %% Band gap according to MPB calculation
set(ha,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
hold on
hb = area([ffield(1)*10^-12 332.5], [5000 5000]); %% Band gap according to MPB calculation
hc = area([(fg*10^(-12))-0.1,(fg*10^(-12))+1],[5000 5000]);
set(hb,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
set(hc,'Facecolor',[1 0.85 1],'Edgecolor','none')
set (gca, 'Yscale', 'log');
set(gca,'xticklabel',{[]}) 
set(gca,'YMinorTick','on')
h1=semilogy(ffield*10^(-12),abs((Jij_yy-Jij0)./Gamma0),'k-','linewidth',1.2)
h2=line([351.7250 351.7250], [0.01 1000])
h3=line([380.6528 380.6528], [0.01 1000])
hold off
xlim([ffield(1)*10^(-12) ffield(end)*10^(-12)])
ylim([0.01 5000])
set(h2,'color',[0.9 0.2 0.2],'linestyle','--','linewidth',0.8)
set(h3,'color',[0.2 0.2 0.9],'linestyle','--','linewidth',1)
text(349.5,2500,'$\nu_{D2}$','Interpreter','latex','Fontsize',18,'Color',[0.9 0.2 0.2])
text(378,2500,'$\nu_{\rm trap}$','Interpreter','latex','Fontsize',18,'Color',[0.2 0.2 0.9])
grid on
set(gca,'XTick',[300:10:400],'Fontsize',12)
% set(gca,'Position',[0.08 0.68 0.6 0.28])
Ylabel=ylabel('$\mid \Delta_{\rm vdW}^{\rm 1D}/\Gamma_{\rm vac} \mid$ ','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.01,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.5,0,0])
set(gca,'Layer','top')

subplot(4,2,4,'position',[0.7 0.53 0.22 0.21])
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),Gamma_yy./Gamma0,'k-','linewidth',1.5)
hold on
hold off
xlim([-0.1 1])
ylim([0.4 10000])
grid on
%set(h2,'color',[0.9 0.3 0.3],'linestyle','--','linewidth',1)
set(gca,'xticklabel',{[]}) 
set(gca,'YMinorTick','on')
text(351.45,500,'$\nu_{D2}$','Interpreter','latex','Fontsize',18,'Color',[0.9 0.3 0.3])
set(gca,'Fontsize',12)
% Xlabel=xlabel('Frequency (THz)','Interpreter','Latex','Fontsize',17)
% Ylabel=ylabel('$\Gamma_{\rm total}/\Gamma_{\rm vac}$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.1,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')

subplot(4,2,3,'Position',[0.08 0.53 0.58 0.21])
ha = area([384.1 ffield(end)*10^-12], [10000 10000]); %% Band gap according to MPB calculation
set(ha,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
hold on
hb = area([ffield(1)*10^-12 332.5], [10000 10000]); %% Band gap according to MPB calculation
hc = area([(fg*10^(-12))-0.1,(fg*10^(-12))+1],[10000 10000]);
set(hb,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
set(hc,'Facecolor',[1 0.85 1],'Edgecolor','none')
set (gca, 'Yscale', 'log');
set(gca,'xticklabel',{[]}) 
set(gca,'YMinorTick','on')
h1=semilogy(ffield*10^(-12),Gamma_yy./Gamma0,'k-','linewidth',1.2)
h2=line([351.7250 351.7250], [0.4 10000])
h3=line([380.6528 380.6528], [0.4 10000])
hold off
xlim([ffield(1)*10^(-12) ffield(end)*10^(-12)])
ylim([0.4 10000])
set(h2,'color',[0.9 0.2 0.2],'linestyle','--','linewidth',1)
set(h3,'color',[0.2 0.2 0.9],'linestyle','--','linewidth',1)
grid on
set(gca,'XTick',[300:10:400],'Fontsize',12)
Ylabel=ylabel('$\Gamma_{\rm total}/\Gamma_{\rm vac}$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.01,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.5,0,0])
set(gca,'Layer','top')

subplot(4,2,6,'position',[0.7 0.30 0.22 0.21])
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma_yy),'k-','linewidth',1.5)
hold on
hold off
xlim([-0.1 1])
ylim([0.01 100])
grid on
set(gca,'Fontsize',12)
% Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
set(gca,'xticklabel',{[]})
set(gca,'YMinorTick','on')
% Ylabel=ylabel('$\mid J_{\rm 1D}/\Gamma_{\rm total} \mid$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')

subplot(4,2,5,'Position',[0.08 0.30 0.58 0.21])
ha = area([384.1 ffield(end)*10^-12], [200 200]); %% Band gap according to MPB calculation
set(ha,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
hold on
hb = area([ffield(1)*10^-12 332.5], [200 200]); %% Band gap according to MPB calculation
hc = area([(fg*10^(-12))-0.1,(fg*10^(-12))+1],[200 200]);
set(hb,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
set(hc,'Facecolor',[1 0.85 1],'Edgecolor','none')
set (gca, 'Yscale', 'log');
h1=semilogy(ffield*10^(-12),abs((Jij_yy-Jij0)./Gamma_yy),'k-','linewidth',1.2)
h2=line([351.7250 351.7250], [0.01 200])
h3=line([380.6528 380.6528], [0.01 200])
hold off
xlim([ffield(1)*10^(-12) ffield(end)*10^(-12)])
ylim([0.01 100])
set(h2,'color',[0.9 0.2 0.2],'linestyle','--','linewidth',1)
set(h3,'color',[0.2 0.2 0.9],'linestyle','--','linewidth',1)
set(gca,'xticklabel',{[]}) 
set(gca,'YMinorTick','on')
grid on
set(gca,'XTick',[300:10:400],'Fontsize',12)
% Xlabel=xlabel('Frequency (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\mid \Delta_{\rm vdW}^{\rm 1D}/\Gamma_{\rm total} \mid$ ','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.5,0,0])
set(gca,'Layer','top')

subplot(4,2,8,'position',[0.7 0.07 0.22 0.21])
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma1D),'k-','linewidth',1.5)
hold on
hold off
xlim([-0.1 1])
ylim([0.01 1000])
grid on
set(gca,'Fontsize',12)
set(gca,'YMinorTick','on')
Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
% Ylabel=ylabel('$\mid J_{\rm 1D}/\Gamma_{\rm total} \mid$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')

subplot(4,2,7,'Position',[0.08 0.07 0.58 0.21])
ha = area([384.1 ffield(end)*10^-12], [1000 1000]); %% Band gap according to MPB calculation
set(ha,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
hold on
hb = area([ffield(1)*10^-12 332.5], [1000 1000]); %% Band gap according to MPB calculation
hc = area([(fg*10^(-12))-0.1,(fg*10^(-12))+1],[1000 1000]);
set(hb,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
set(hc,'Facecolor',[1 0.85 1],'Edgecolor','none')
set (gca, 'Yscale', 'log');
h1=semilogy(ffield*10^(-12),abs((Jij_yy-Jij0)./Gamma1D),'k-','linewidth',1.2)
h2=line([351.7250 351.7250], [0.01 1000])
h3=line([380.6528 380.6528], [0.01 1000])
hold off
xlim([ffield(1)*10^(-12) ffield(end)*10^(-12)])
ylim([0.01 1000])
set(h2,'color',[0.9 0.2 0.2],'linestyle','--','linewidth',1)
set(h3,'color',[0.2 0.2 0.9],'linestyle','--','linewidth',1)
grid on
set(gca,'XTick',[300:10:400],'Fontsize',12)
set(gca,'YMinorTick','on')
Xlabel=xlabel('Frequency (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\mid \Delta_{\rm vdW}^{\rm 1D}/\Gamma_{\rm 1D} \mid$ ','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.5,0,0])
set(gca,'Layer','top')
saveas(A,['J1DGamma.png'])

%% Zoom Gamma

B=figure(2)
set(B,'Position',[1000 800 600 500]);
hb = area([fg*10^(-12)-0.3 fg*10^(-12)+0.4], [10000 10000]); %% Band gap according to MPB calculation
set(hb,'Facecolor',[0.95 0.95 0.8],'Edgecolor','none')
set (gca, 'Yscale', 'log');
hold on
h=semilogy(ffield*10^(-12),Gamma_yy./Gamma0,'k-','linewidth',1.5)
h2=line([351.7250 351.7250], [1 10000])
hold off
xlim([350 352])
ylim([1 10^4])
grid on
set(h2,'color',[0.9 0.3 0.3],'linestyle','--','linewidth',1)
text(351.45,500,'$\nu_{D2}$','Interpreter','latex','Fontsize',18,'Color',[0.9 0.3 0.3])
set(gca,'XTick',[348:1:353],'Fontsize',12)
Xlabel=xlabel('Frequency (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\Gamma_{\rm total}/\Gamma_{\rm vac}$','Interpreter','Latex','Fontsize',18)
set(gca,'Position',[0.12 0.15 0.85 0.8])
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.1,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')
saveas(B,['Green_yy.png'])

F=figure(10)
set(F,'Position',[1000 800 600 500]);
h=plot(ffield*10^(-12),Gamma_yy./Gamma0,'k-','linewidth',1.5)
hold on
xlim([350.50 350.7])
% ylim([0 700])
grid on
set(gca,'XTick',[350.5:0.05:350.7],'Fontsize',12)
Xlabel=xlabel('Frequency (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\Gamma_{\rm total}/\Gamma_{\rm vac}$','Interpreter','Latex','Fontsize',18)
set(gca,'Position',[0.12 0.15 0.85 0.8])
[M1,I1]=min(abs(ffield*(10^-12)-350.5));
[M2,I2]=min(abs(ffield*(10^-12)-350.7));
func=@(c,x)c(1)./(((x-c(2)).^2)+(c(3)).^2);
c0=[1 350.6 0.001];
bn=lsqcurvefit(func,c0,ffield(I1:I2)*10^(-12),Gamma_yy(I1:I2)./Gamma0(I1:I2))
W=linspace(ffield(I1)*10^(-12),ffield(I2)*10^(-12),1000);
h1=plot(W,func(bn,W),'r--','linewidth',2);
hold off
title('Lorentzian fitting')
text(350.61,6500,'$f(\nu)=\frac{a_{0}}{((\nu-a_{1})^2-a_{2}^2)}$','Interpreter','latex','fontsize',16)
text(350.61,6000,['$a_{0}=$ ' num2str(bn(1)) ],'Interpreter','latex','fontsize',13)
text(350.61,5500,['$a_{1}=$ ' num2str(bn(2)) ],'Interpreter','latex','fontsize',13)
text(350.61,5000,['$a_{2}=$ ' num2str(abs(bn(3))) ],'Interpreter','latex','fontsize',13)
text(350.61,4500,['linewidth = ' num2str(abs(2*bn(3))) 'THz'],'Interpreter','latex','fontsize',13)
text(350.61,4000,['$Q=\nu_{0}/\delta \nu=$' num2str(bn(2)/abs(2*bn(3)))],'Interpreter','latex','fontsize',13)
saveas(F,['LorentzianFit.png'])

%% Zoom Jij/Gamma_total

C=figure(3)
set(C,'Position',[1000 800 600 500]);
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma_yy),'k-','linewidth',1.5)
xlim([-0.2 0.3])
ylim([0.01 100])
grid on
set(gca,'Fontsize',12)
set(gca,'Position',[0.12 0.15 0.85 0.8])
Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\mid \Delta_{\rm vdW}^{\rm 1D}/\Gamma_{\rm total} \mid$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')
saveas(C,['J1D_GammaTotal.png'])

%% Zoom Jij/Gamma1D

D=figure(4)
set(D,'Position',[1000 800 600 500]);
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma1D),'k-','linewidth',1.2)
xlim([-0.2 0.3])
grid on
set(gca,'Fontsize',12)
set(gca,'Position',[0.12 0.15 0.85 0.8])
Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\mid \Delta_{\rm vdW}^{\rm 1D}/\Gamma_{\rm 1D} \mid$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')
saveas(D,['J1D_Gamma1D.png'])

%% Nonlocal Green's function
YareaM=repmat(-1,length(Xplot),1);
YareaP=repmat(1,length(Xplot),1);
x2=[Xplot'*10^6;flipud(Xplot'*10^6)];
inBetween=[YareaM;flipud(YareaP)];
Index=Index+200

E=figure(5)
set(E,'PaperPositionMode','auto')
subplot(6,1,1,'Position',[0.13 0.85 0.8 0.13])
set(E,'Position',[1000 1000 1000 1000]);
ha = fill(x2,inBetween,[0.85 0.85 0.85]); 
hold on
h=plot(Xplot*10^6,(NLoJij(:,1)-NLoJij0(:,1))/NLoGamma(50001,1),'k-');
h=plot(Xplot(Index)*10^6,(NLoJij(Index,1)-NLoJij0(Index,1))/NLoGamma(50001,1),'b.');
set(ha,'Edgecolor',[0.85 0.85 0.85])
xlim([Xplot(1)*10^6,Xplot(end)*10^6])
%ylim([-0.1 0.1])
set(gca,'XTick',[-16:4:16])
set(gca,'Fontsize',12)
ylabel([])
set(gca,'xticklabel',{[]})
text(10,0.065,'$\Delta_{\rm e}=0$ THz','Interpreter','Latex','Fontsize',16')
set(gca,'Layer','top')
grid on

subplot(6,1,2,'Position',[0.13 0.7 0.8 0.13])
ha = fill(x2,inBetween,[0.85 0.85 0.85]); 
hold on
h=plot(Xplot*10^6,(NLoJij(:,2)-NLoJij0(:,2))/NLoGamma(50001,2),'k-');
h=plot(Xplot(Index)*10^6,(NLoJij(Index,2)-NLoJij0(Index,2))/NLoGamma(50001,2),'b.');
set(ha,'Edgecolor','none')
xlim([Xplot(1)*10^6,Xplot(end)*10^6])
set(gca,'XTick',[-16:4:16])
set(gca,'Fontsize',12)
ylabel([])
set(gca,'xticklabel',{[]}) 
text(10,13,'$\Delta_{\rm e}=0.1$ THz','Interpreter','Latex','Fontsize',16')
set(gca,'Layer','top')
grid on

subplot(6,1,3,'Position',[0.13 0.55 0.8 0.13])
ha = fill(x2,inBetween,[0.85 0.85 0.85]); 
hold on
h1=plot(Xplot*10^6,(NLoJij(:,3)-NLoJij0(:,3))/NLoGamma(50001,3),'k-');
h=plot(Xplot(Index)*10^6,(NLoJij(Index,3)-NLoJij0(Index,3))/NLoGamma(50001,3),'b.');
xlim([Xplot(1)*10^6,Xplot(end)*10^6])
set(gca,'XTick',[-16:4:16])
set(gca,'Fontsize',12)
ylabel('$\Delta_{\rm vdW}^{\rm ij}/\Gamma_{\rm total}$','Interpreter','Latex','Fontsize',16')
set(ha,'Edgecolor','none')
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,10,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,11.5,0])
set(gca,'Layer','top')
set(gca,'xticklabel',{[]}) 
text(10,6.5,'$\Delta_{\rm e}=0.2$ THz','Interpreter','Latex','Fontsize',16')
grid on

subplot(6,1,4,'Position',[0.13 0.4 0.8 0.13])
ha = fill(x2,inBetween,[0.85 0.85 0.85]); 
hold on
h2=plot(Xplot*10^6,(NLoJij(:,4)-NLoJij0(:,4))/NLoGamma(50001,4),'k-');
h=plot(Xplot(Index)*10^6,(NLoJij(Index,4)-NLoJij0(Index,4))/NLoGamma(50001,4),'b.');
xlim([Xplot(1)*10^6,Xplot(end)*10^6])
set(ha,'Edgecolor','none')
set(gca,'XTick',[-16:4:16])
set(gca,'xticklabel',{[]}) 
set(gca,'Fontsize',12)
ylabel([])
text(10,6.5,'$\Delta_{\rm e}=0.3$ THz','Interpreter','Latex','Fontsize',16')
set(gca,'Layer','top')
grid on

subplot(6,1,5,'Position',[0.13 0.25 0.8 0.13])
ha = fill(x2,inBetween,[0.85 0.85 0.85]); 
hold on
h3=plot(Xplot*10^6,(NLoJij(:,5)-NLoJij0(:,5))/NLoGamma(50001,5),'k-');
h=plot(Xplot(Index)*10^6,(NLoJij(Index,5)-NLoJij0(Index,5))/NLoGamma(50001,5),'b.');
set (gca, 'Yscale', 'linear');
set(gca,'Fontsize',12)
xlim([Xplot(1)*10^6,Xplot(end)*10^6])
set(ha,'Edgecolor','none')
set(gca,'xticklabel',{[]}) 
xlabel([])
ylabel([])
text(10,6.5,'$\Delta_{\rm e}=0.4$ THz','Interpreter','Latex','Fontsize',16')
set(gca,'Layer','top')
grid on

subplot(6,1,6,'Position',[0.13 0.1 0.8 0.13])
ha = fill(x2,inBetween,[0.85 0.85 0.85]); 
hold on
h3=plot(Xplot*10^6,(NLoJij(:,6)-NLoJij0(:,6))/NLoGamma(50001,6),'k-');
h=plot(Xplot(Index)*10^6,(NLoJij(Index,6)-NLoJij0(Index,6))/NLoGamma(50001,6),'b.');
set (gca, 'Yscale', 'linear');
set(gca,'Fontsize',12)
xlim([Xplot(1)*10^6,Xplot(end)*10^6])
set(ha,'Edgecolor','none')
set(gca,'XTick',[-16:4:16],'Fontsize',12)
xlabel('x ($\mu$m)','Interpreter','Latex','Fontsize',16')
ylabel([])
text(10,3.5,'$\Delta_{\rm e}=0.5$ THz','Interpreter','Latex','Fontsize',16')
set(gca,'Layer','top')
grid on
saveas(E,['NLoGF.png'])