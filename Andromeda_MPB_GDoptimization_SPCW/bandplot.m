function []=bandplot(a)
%%%%% constants // Unit conversion %%%%%
c=299792458;
%a=374*10^-9;     % Lattice constant in m
ep0=8.85e-12;
hbar=1.05457148e-34;
h=6.626e-34;
MPB_unit=1;
d2_freq=351.725; % THz
bluemagic_freq=c/(793.5*10^-9)*10^(-12)
SIunit_omega=(c/a)*2*pi*MPB_unit;  %%% frequency in THz
SIunit_energy=SIunit_omega*hbar;
SIunit_prop=(2*pi/a);
SIunit_mass=hbar/(a^2*2*pi*c/a);
SIunit_k0=0.5*2*pi/a;

%% load MAT file %%
load('data.mat');

%% Data preparation %%
data=data{1};
lightcone=data.vectors(:,1).*SIunit_omega./(2*pi).*10^(-12);
kpoints=data.vectors(:,1);
gbands20=data.bands(:,20).*SIunit_omega./(2*pi).*10^(-12);
gbands21=data.bands(:,21).*SIunit_omega./(2*pi).*10^(-12);
gbands22=data.bands(:,22).*SIunit_omega./(2*pi).*10^(-12);
gbands23=data.bands(:,23).*SIunit_omega./(2*pi).*10^(-12);
lowercontinum=data.bands(:,19).*SIunit_omega./(2*pi).*10^(-12);
uppercontinum=data.bands(:,24).*SIunit_omega./(2*pi).*10^(-12);
upperline=repmat(0.6,size(kpoints,1),1).*SIunit_omega./(2*pi).*10^(-12);
D2line=repmat(d2_freq,1,length(kpoints));
Bmagic=repmat(bluemagic_freq,1,length(kpoints));

%% PLOT %%
A=figure(1)
h1=plot(kpoints,gbands20)
hold on
h2=plot(kpoints,gbands21)
h3=plot(kpoints,gbands22)
h5=area(kpoints,lowercontinum)
xx=[kpoints',fliplr(kpoints')];
inBetween1=[lightcone',fliplr(upperline')];
inBetween2=[gbands23',fliplr(upperline')];
h9=fill(xx,inBetween2,[0.8 0.8 0.8],'EdgeColor','none')
h10=fill(xx,inBetween1,[1 1 1],'EdgeColor','none')
h7=plot(kpoints,D2line)
h8=plot(kpoints,Bmagic)
h6=plot(kpoints,lightcone)
xlim([0.42 0.5])
ylim([310 400])
hold off

set(h7,'Color',[0 0 0],'LineStyle','--')
set(h8,'Color',[0 0 0],'LineStyle','--')
set(h1,'Color',[1 0.2 0.1],'Linewidth',2)
% set(h1,'Color',[0.8 0.2 0.4],'Linewidth',1.5)
set(h2,'Color',[0.6 0.6 0.6],'Linewidth',1.5)
set(h3,'Color',[0.1 0.4 1],'Linewidth',2)
% set(h3,'Color',[0.2 0.2 0.7],'Linewidth',1.5)
set(h6,'Color',[0 0 0],'Linewidth',1.5)
set(h5,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
set(gca,'Layer','top')
set(gca,'ytick',[320 340 360 380 400])
set(gca,'Fontsize',12)
xlabel(['Propagation Constant ($$k_{x}a_{0}/2\pi$$)'],'Interpreter','Latex','FontSize',16)
ylabel(['Frequency (THz)'],'Interpreter','Latex','FontSize',16)
hText = text(0.425, 356,'$$\nu_{D2}$$','Interpreter','Latex','FontSize',16);
set(gcf, 'PaperSize', [4 2]);
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,1,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.0005,0,0])
set(gca,'Position',[0.15 0.15 0.7 0.8])
saveas(A,'Band.png');

B=figure(2)
h1=plot(kpoints,data.bands(:,20))
hold on
h2=plot(kpoints,data.bands(:,21))
h3=plot(kpoints,data.bands(:,22))
h5=area(kpoints,data.bands(:,19))
xx=[kpoints',fliplr(kpoints')];
inBetween1=[data.vectors(:,1)',fliplr(repmat(0.6,size(kpoints,1),1)')];
inBetween2=[data.bands(:,23)',fliplr(repmat(0.6,size(kpoints,1),1)')];
h9=fill(xx,inBetween2,[0.8 0.8 0.8],'EdgeColor','none')
h10=fill(xx,inBetween1,[1 1 1],'EdgeColor','none')
h7=plot(kpoints,D2line/(SIunit_omega./(2*pi).*10^(-12)))
h8=plot(kpoints,Bmagic/(SIunit_omega./(2*pi).*10^(-12)))
h6=plot(kpoints,data.vectors(:,1))
xlim([0.42 0.5])
ylim([0.39 0.5])
hold off

set(h7,'Color',[0 0 0],'LineStyle','--')
set(h8,'Color',[0 0 0],'LineStyle','--')
set(h1,'Color',[1 0.2 0.1],'Linewidth',2)
set(h2,'Color',[0.6 0.6 0.6],'Linewidth',1.5)
set(h3,'Color',[0.1 0.4 1],'Linewidth',2)
set(h6,'Color',[0 0 0],'Linewidth',1.5)
set(h5,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
set(gca,'Layer','top')
set(gca,'Fontsize',12)
xlabel(['Propagation Constant ($$k_{x}a_{0}/2\pi$$)'],'Interpreter','Latex','FontSize',16)
ylabel(['Normalized Frequency $$(\omega a_{0}/2 \pi c)$$'],'Interpreter','Latex','FontSize',16)
hText = text(0.425, 0.4426,'$$\nu_{D2}$$','Interpreter','Latex','FontSize',16);
set(gcf, 'PaperSize', [4 2]);
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.0005,0,0])
set(gca,'Position',[0.15 0.15 0.7 0.8])
saveas(B,'NormBand.png');


C=figure(3)
h1=plot(kpoints,data.bands(:,20))
hold on
h2=plot(kpoints,data.bands(:,21))
h3=plot(kpoints,data.bands(:,22))
h5=area(kpoints,data.bands(:,19))
xx=[kpoints',fliplr(kpoints')];
inBetween1=[data.vectors(:,1)',fliplr(repmat(0.6,size(kpoints,1),1)')];
inBetween2=[data.bands(:,23)',fliplr(repmat(0.6,size(kpoints,1),1)')];
h9=fill(xx,inBetween2,[0.8 0.8 0.8],'EdgeColor','none')
h10=fill(xx,inBetween1,[1 1 1],'EdgeColor','none')
h7=plot(kpoints,D2line/(SIunit_omega./(2*pi).*10^(-12)))
h8=plot(kpoints,Bmagic/(SIunit_omega./(2*pi).*10^(-12)))
h6=plot(kpoints,data.vectors(:,1))
xlim([0.42 0.5])
ylim([0.42 0.44])
hold off

set(h7,'Color',[0 0 0],'LineStyle','--')
set(h8,'Color',[0 0 0],'LineStyle','--')
set(h1,'Color',[1 0.2 0.1],'Linewidth',2)
set(h2,'Color',[0.6 0.6 0.6],'Linewidth',1.5)
set(h3,'Color',[0.1 0.4 1],'Linewidth',2)
set(h6,'Color',[0 0 0],'Linewidth',1.5)
set(h5,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
set(gca,'Layer','top')
set(gca,'Fontsize',12)
xlabel(['Propagation Constant ($$k_{x}a_{0}/2\pi$$)'],'Interpreter','Latex','FontSize',16)
ylabel(['Normalized Frequency $$(\omega a_{0}/2 \pi c)$$'],'Interpreter','Latex','FontSize',16)
hText = text(0.425, 0.4426,'$$\nu_{D2}$$','Interpreter','Latex','FontSize',16);
set(gcf, 'PaperSize', [4 2]);
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.0005,0,0])
set(gca,'Position',[0.15 0.15 0.7 0.8])
saveas(C,'ZoomNormband.png')
end
