%% contants
clear all
close all

clight=299792458;
a=374*10^-9;     % Lattice constant in m
eps0=8.85e-12;
hbar=1.05457148e-34;
h=6.626e-34;

%% Data load
load('GFresulty.mat')
ffieldy=GFresult.f;
Efieldy=GFresult.E;
X=GFresult.x;
load('GFresult0y.mat')
Efield0y=GFresult0.E;
ffield0y=GFresult0.f;
X0=GFresult0.x;

load('muresulty.mat')
load('GreenNf2D.mat')
load('GreenN2D.mat')
load('Green0f2D.mat')
load('Green02D.mat')

Xplot=linspace(X(1),X(end),100000);
ffield=linspace(ffieldy(1),ffieldy(end),100000);
InG0f=interp1(ffield0y,squeeze(G0f(2,2,:)),ffield,'linear');
InGNfyy=interp1(ffieldy,squeeze(GNf(2,2,:)),ffield,'linear');
InG0=interp1(X0,squeeze(G0(2,2,:)),Xplot,'linear');
InGNyy=interp1(X,squeeze(GN(2,2,:)),Xplot,'linear');

[M I]=max(squeeze(imag(GNf(2,2,:))));
fg=ffieldy(I);
delta=(ffieldy-fg);
deltav=(ffield0y-fg);
[m0 q0]=min(abs(delta));
[mv0 qv0]=min(abs(deltav));

delta1=0.1*10^12; [m1 q1]=min(abs(delta-delta1)); [mv1 qv1]=min(abs(deltav-delta1));
delta2=0.2*10^12; [m2 q2]=min(abs(delta-delta2)); [mv2 qv2]=min(abs(deltav-delta2));
delta3=0.3*10^12; [m3 q3]=min(abs(delta-delta3)); [mv3 qv3]=min(abs(deltav-delta3));
delta4=0.4*10^12; [m4 q4]=min(abs(delta-delta4)); [mv4 qv4]=min(abs(deltav-delta4));
delta5=0.5*10^12; [m5 q5]=min(abs(delta-delta5)); [mv5 qv5]=min(abs(deltav-delta5));
delta6=0.6*10^12; [m6 q6]=min(abs(delta-delta6)); [mv6 qv6]=min(abs(deltav-delta6));

w0=2*pi*ffieldy(q0); wv0=2*pi*ffield0y(qv0);
w1=2*pi*ffieldy(q1); wv1=2*pi*ffield0y(qv1);
w2=2*pi*ffieldy(q2); wv2=2*pi*ffield0y(qv2);
w3=2*pi*ffieldy(q3); wv3=2*pi*ffield0y(qv3);
w4=2*pi*ffieldy(q4); wv4=2*pi*ffield0y(qv4);
w5=2*pi*ffieldy(q5); wv5=2*pi*ffield0y(qv5);
w6=2*pi*ffieldy(q6); wv6=2*pi*ffield0y(qv6);

irun=0;
q=[q0 q1 q2 q3 q4 q5 q6];
qv=[qv0 qv1 qv2 qv3 qv4 qv5 qv6];
w=[w0 w1 w2 w3 w4 w5 w6];
wv=[wv0 wv1 wv2 wv3 wv4 wv5 wv6];

for j=1:length(q)
    Gyy(:,j)=squeeze((Efieldy(:,2,q(j))))*clight^2*eps0/w(j)^2/mu;
    Gyy0(:,j)=squeeze((Efield0y(:,2,qv(j))))*clight^2*eps0/wv(j)^2/mu;
end

for i=1:length(q)
    InGyy(:,i)=interp1(X,squeeze(Gyy(:,i)),Xplot);
    InGyy0(:,i)=interp1(X0,squeeze(Gyy0(:,i)),Xplot);
end

%% Constant
clight=299792458;
mu0=4*pi*10^(-7);
eps0=8.85e-12;
hbar=1.05457148e-34;
h=2*pi*hbar;
fcs_d1=335.116*10^12;
fcs_d2=351.725*10^12; % Units in Hz
omega_d1=2*pi*fcs_d1;
omega_d2=2*pi*fcs_d2;
ec=1.6022*10^-19;
a0=5.291*10^-11;
JdJD1=4.489*ec*a0; % Jonathan's thesis (A.19)
JdJD2=6.324*ec*a0; % Jonathan's thesis (A.20) : already includes jj_factor
jj_factor=0.5; % D2 transition

%% Gamma & Jij
Gamma_yy=(JdJD2*jj_factor)^2*(2*omega_d2^2/(hbar*eps0*clight^2))*(imag(InGNfyy));
Jij_yy=(JdJD2*jj_factor)^2*(omega_d2^2/(hbar*eps0*clight^2))*(real(InGNfyy));
Gamma0=(JdJD2*jj_factor)^2*(2*omega_d2^2/(hbar*eps0*clight^2))*(imag(InG0f));
Jij0=(JdJD2*jj_factor)^2*(omega_d2^2/(hbar*eps0*clight^2))*(real(InG0f));
Gamma1D=Gamma_yy-Gamma_yy(34000);

NLoGamma=(JdJD2*jj_factor)^2*(2*omega_d2^2/(hbar*eps0*clight^2))*(imag(InGyy));
NLoJij=(JdJD2*jj_factor)^2*(omega_d2^2/(hbar*eps0*clight^2))*(real(InGyy));
NLoGamma0=(JdJD2*jj_factor)^2*(2*omega_d2^2/(hbar*eps0*clight^2))*(imag(InGyy0));
NLoJij0=(JdJD2*jj_factor)^2*(omega_d2^2/(hbar*eps0*clight^2))*(real(InGyy0));

save('Gamma_tot.mat','Gamma_yy')
save('Gamma_vac.mat','Gamma0')
save('Jij.mat','Jij_yy')
save('Jij_vac.mat','Jij0')
save('Gamma1D.mat','Gamma1D')
save('freq.mat','ffield')
save('NLoGamma.mat','NLoGamma')
save('NLoGamma0.mat','NLoGamma0')
save('NLoJij.mat','NLoJij')
save('NLoJij0.mat','NLoJij0')
save('Xplot.mat','Xplot')
%% plot
A=figure(1)
set(A,'PaperPositionMode','auto')
set(A,'Position',[1000 800 1400 1000]);
H=subplot(3,2,2,'Position',[0.7 0.68 0.25 0.28])
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma0),'k-','linewidth',1.5)
xlim([-0.3 0.5])
ylim([3 1000])
grid on
set(gca,'Fontsize',12)
set(gca,'xticklabel',{[]}) 
% Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
% Ylabel=ylabel('$\mid J_{\rm 1D} \mid$ (MHz) ','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')

H=subplot(3,2,1,'Position',[0.08 0.68 0.58 0.28])
ha = area([386.1287 ffield(end)*10^-12], [1000 1000]); %% Band gap according to MPB calculation
set(ha,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
hold on
hb = area([ffield(1)*10^-12 334.3634], [1000 1000]); %% Band gap according to MPB calculation
hc = area([(fg*10^(-12))-0.3,(fg*10^(-12))+0.5],[1000 1000]);
set(hb,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
set(hc,'Facecolor',[1 0.85 1],'Edgecolor','none')
set (gca, 'Yscale', 'log');
set(gca,'xticklabel',{[]}) 
h1=semilogy(ffield*10^(-12),abs((Jij_yy-Jij0)./Gamma0),'k-','linewidth',1.2)
h2=line([351.7250 351.7250], [0.01 1000])
h3=line([380.6528 380.6528], [0.01 1000])
hold off
xlim([ffield(1)*10^(-12) ffield(end)*10^(-12)])
ylim([0.01 1000])
set(h2,'color',[0.9 0.2 0.2],'linestyle','--','linewidth',0.8)
set(h3,'color',[0.2 0.2 0.9],'linestyle','--','linewidth',1)
text(349.5,2000,'$\nu_{D2}$','Interpreter','latex','Fontsize',18,'Color',[0.9 0.2 0.2])
text(378,2000,'$\nu_{\rm trap}$','Interpreter','latex','Fontsize',18,'Color',[0.2 0.2 0.9])
grid on
set(gca,'XTick',[300:10:400],'Fontsize',12)
% set(gca,'Position',[0.08 0.68 0.6 0.28])
Ylabel=ylabel('$\mid J_{\rm 1D}/\Gamma_{\rm vac} \mid$ ','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.01,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.5,0,0])
set(gca,'Layer','top')

subplot(3,2,4,'position',[0.7 0.38 0.25 0.28])
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),Gamma_yy./Gamma0,'k-','linewidth',1.5)
h2=line([351.7250 351.7250], [1 2000])
hold off
xlim([-0.3 0.5])
ylim([0.4 2000])
grid on
set(h2,'color',[0.9 0.3 0.3],'linestyle','--','linewidth',1)
set(gca,'xticklabel',{[]}) 
text(351.45,500,'$\nu_{D2}$','Interpreter','latex','Fontsize',18,'Color',[0.9 0.3 0.3])
set(gca,'Fontsize',12)
% Xlabel=xlabel('Frequency (THz)','Interpreter','Latex','Fontsize',17)
% Ylabel=ylabel('$\Gamma_{\rm total}/\Gamma_{\rm vac}$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.1,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')

subplot(3,2,3,'Position',[0.08 0.38 0.58 0.28])
ha = area([386.1287 ffield(end)*10^-12], [10000 10000]); %% Band gap according to MPB calculation
set(ha,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
hold on
hb = area([ffield(1)*10^-12 334.3634], [10000 10000]); %% Band gap according to MPB calculation
hc = area([(fg*10^(-12))-0.3,(fg*10^(-12))+0.5],[10000 10000]);
set(hb,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
set(hc,'Facecolor',[1 0.85 1],'Edgecolor','none')
set (gca, 'Yscale', 'log');
set(gca,'xticklabel',{[]}) 
h1=semilogy(ffield*10^(-12),Gamma_yy./Gamma0,'k-','linewidth',1.2)
h2=line([351.7250 351.7250], [0.4 10000000])
h3=line([380.6528 380.6528], [0.4 10000000])
hold off
xlim([ffield(1)*10^(-12) ffield(end)*10^(-12)])
ylim([0.4 2000])
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
saveas(A,['Gamma.png'])

subplot(3,2,6,'position',[0.7 0.08 0.25 0.28])
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma_yy),'k-','linewidth',1.5)
hold on
h3=line([0.1, 0.1], [0.01 100],'color',[1 0 0],'linestyle','--')
h4=line([0.2, 0.2], [0.01 100],'color',[1 0 1],'linestyle','--')
h5=line([0.3 0.3], [0.01 100],'color',[0 1 0],'linestyle','--')
h6=line([0.4 0.4], [0.01 100],'color',[0 0 1],'linestyle','--')
hold off
text(0.11,0.02,'$\nu_{1}$','Interpreter','latex','Fontsize',18,'Color',[1 0 0])
text(0.21,0.02,'$\nu_{2}$','Interpreter','latex','Fontsize',18,'Color',[1 0 1])
text(0.31,0.02,'$\nu_{3}$','Interpreter','latex','Fontsize',18,'Color',[0 1 0])
text(0.41,0.02,'$\nu_{4}$','Interpreter','latex','Fontsize',18,'Color',[0 0 1])
xlim([-0.3 0.5])
ylim([0.01 100])
grid on
set(gca,'Fontsize',12)
Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
% Ylabel=ylabel('$\mid J_{\rm 1D}/\Gamma_{\rm total} \mid$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')

subplot(3,2,5,'Position',[0.08 0.08 0.58 0.28])
ha = area([386.1287 ffield(end)*10^-12], [200 200]); %% Band gap according to MPB calculation
set(ha,'Facecolor',[0.85 0.85 0.85],'Edgecolor','none')
hold on
hb = area([ffield(1)*10^-12 334.3634], [200 200]); %% Band gap according to MPB calculation
hc = area([(fg*10^(-12))-0.3,(fg*10^(-12))+0.5],[200 200]);
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
grid on
set(gca,'XTick',[300:10:400],'Fontsize',12)
Xlabel=xlabel('Frequency (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\mid J_{\rm 1D}/\Gamma_{\rm total} \mid$ ','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.5,0,0])
set(gca,'Layer','top')
saveas(A,['J1DGamma.png'])


%%
B=figure(2)
set(B,'Position',[1000 800 600 500]);
hb = area([fg*10^(-12)-0.2 fg*10^(-12)+0.3], [1000 1000]); %% Band gap according to MPB calculation
set(hb,'Facecolor',[0.95 0.95 0.8],'Edgecolor','none')
set (gca, 'Yscale', 'log');
hold on
h=semilogy(ffield*10^(-12),Gamma_yy./Gamma0,'k-','linewidth',1.5)
h2=line([351.7250 351.7250], [1 1000])
hold off
xlim([348 351])
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


%%
C=figure(3)
set(C,'Position',[1000 800 600 500]);
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma_yy),'k-','linewidth',1.5)
hold on
h3=line([0.05, 0.05], [0.01 10],'color',[1 0 0],'linestyle','--')
h4=line([0.1, 0.1], [0.01 10],'color',[1 0 1],'linestyle','--')
h5=line([0.15 0.15], [0.01 10],'color',[0 1 0],'linestyle','--')
h6=line([0.2 0.2], [0.01 10],'color',[0 0 1],'linestyle','--')
hold off
text(0.06,0.02,'$\nu_{1}$','Interpreter','latex','Fontsize',18,'Color',[1 0 0])
text(0.11,0.02,'$\nu_{2}$','Interpreter','latex','Fontsize',18,'Color',[1 0 1])
text(0.16,0.02,'$\nu_{3}$','Interpreter','latex','Fontsize',18,'Color',[0 1 0])
text(0.21,0.02,'$\nu_{4}$','Interpreter','latex','Fontsize',18,'Color',[0 0 1])
xlim([-0.2 0.3])
grid on
set(gca,'Fontsize',12)
set(gca,'Position',[0.12 0.15 0.85 0.8])
Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\mid J_{\rm 1D}/\Gamma_{\rm total} \mid$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')
saveas(C,['J1D_GammaTotal.png'])



%%
D=figure(4)
set(D,'Position',[1000 800 600 500]);
set (gca, 'Yscale', 'log');
h=semilogy((ffield-fg)*10^(-12),abs((Jij_yy-Jij0)./Gamma1D),'k-','linewidth',1.2)
hold on
h3=line([0.05, 0.05], [0.01 100],'color',[1 0 0],'linestyle','--')
h4=line([0.1, 0.1], [0.01 100],'color',[1 0 1],'linestyle','--')
h5=line([0.15 0.15], [0.01 100],'color',[0 1 0],'linestyle','--')
h6=line([0.2 0.2], [0.01 100],'color',[0 0 1],'linestyle','--')
hold off
text(0.06,0.02,'$\nu_{1}$','Interpreter','latex','Fontsize',18,'Color',[1 0 0])
text(0.11,0.02,'$\nu_{2}$','Interpreter','latex','Fontsize',18,'Color',[1 0 1])
text(0.16,0.02,'$\nu_{3}$','Interpreter','latex','Fontsize',18,'Color',[0 1 0])
text(0.21,0.02,'$\nu_{4}$','Interpreter','latex','Fontsize',18,'Color',[0 0 1])
xlim([-0.2 0.3])
grid on
set(gca,'Fontsize',12)
set(gca,'Position',[0.12 0.15 0.85 0.8])
Xlabel=xlabel('$\Delta$ (THz)','Interpreter','Latex','Fontsize',17)
Ylabel=ylabel('$\mid J_{\rm 1D}/\Gamma_{\rm 1D} \mid$','Interpreter','Latex','Fontsize',18)
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0,0,0])
set(gca,'Layer','top')
saveas(D,['J1D_Gamma1D.png'])

%%
E=figure(5)
set(E,'Position',[1000 800 600 500]);
set (gca, 'Yscale', 'linear');
h=plot(Xplot*10^6,(NLoJij(:,1)-NLoJij0(:,1))/NLoGamma(25001,2),'r-');
hold on
h1=plot(Xplot*10^6,(NLoJij(:,2)-NLoJij0(:,2))/NLoGamma(25001,3),'b-');
h2=plot(Xplot*10^6,(NLoJij(:,3)-NLoJij0(:,3))/NLoGamma(25001,4),'k-');
h2=plot(Xplot*10^6,(NLoJij(:,4)-NLoJij0(:,4))/NLoGamma(25001,5),'m-');
h2=plot(Xplot*10^6,(NLoJij(:,5)-NLoJij0(:,5))/NLoGamma(25001,6),'c-');
h2=plot(Xplot*10^6,(NLoJij(:,6)-NLoJij0(:,6))/NLoGamma(25001,7),'g-');
hold off
set(gca,'Fontsize',12)
xlim([-8,8])
set(gca,'XTick',[-8:2:8],'Fontsize',12)
xlabel('x ($\mu$ m)','Interpreter','Latex','Fontsize',16')
ylabel('$J_{1 \rm D}$','Interpreter','Latex','Fontsize',16')
% legend('0.1 THz','0.5 THz','1 THz','2 THz', '4 THz', '10 THz')
legend('2 THz', '15 THz', '20 THz')

grid on
saveas(E,['NLoGreen.png'])

F=figure(6)
set(F,'Position',[1000 800 600 500]);
set (gca,'Yscale','linear');
h=plot(Xplot*10^6,NLoGamma(:,2),'r-');
hold on
h1=plot(Xplot*10^6,NLoGamma(:,3),'b-');
h2=plot(Xplot*10^6,NLoGamma(:,4),'k-');
hold off
set(gca,'Fontsize',12)
set(gca,'XTick',[-8:2:8],'Fontsize',12)
xlabel('x ($\mu$ m)','Interpreter','Latex','Fontsize',16')
ylabel('$\Gamma_{total}$','Interpreter','Latex','Fontsize',16')

