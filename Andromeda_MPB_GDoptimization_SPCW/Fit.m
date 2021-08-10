%%%%% constants // Unit conversion %%%%%
c=299792458;
a=374*10^-9;     % Lattice constant in m
ep0=8.85e-12;
hbar=1.05457148e-34;
h=6.626e-34;
d2_freq=351.725; % THz
bluemagic_freq=c/(793.5*10^-9)*10^(-12)
SIunit_omega=(c/a)*2*pi;  %%% frequency in THz
SIunit_energy=SIunit_omega*hbar;
SIunit_prop=(2*pi/a);
SIunit_mass=hbar/(a^2*2*pi*c/a);
SIunit_k0=0.5*2*pi/a;
a=373*10^(-9);

%%%% load MAT file %%%%
load('data.mat');

%% Data preparation %%
data=data{1};
lightcone=data.vectors(:,1);
kpoints=data.vectors(:,1);
gbands20=data.bands(:,20);
gbands21=data.bands(:,21);
gbands22=data.bands(:,22);
gbands23=data.bands(:,23);
lowercontinum=data.bands(:,19);
uppercontinum=data.bands(:,24);
upperline=repmat(0.6,size(kpoints,1),1);
D2line=repmat(d2_freq,1,length(kpoints));
Bmagic=repmat(bluemagic_freq,1,length(kpoints));
Group_index=data.xvelocity(end-15:end-1,20).^(-1);

%% PLOT

figure(1)
h=plot(kpoints(end-15:end-1),Group_index,'bo','linewidth',2)
xlabel('Propagation constant ($k_{x}a_{0}/2\pi$)','Interpreter','Latex','FontSize',16)
ylabel('$1/v_{g}$','Interpreter','Latex','FontSize',16)
hold on
dfunc=@(w,x)(2.*w(1).*(x-0.5)+4.*w(2).*(x-0.5).^(3)+6.*w(3).*(x-0.5).^5).^(-1);
w0=[-0.1 -0.01 -0.01];
an=lsqcurvefit(dfunc,w0,kpoints(end-15:end-1),Group_index(1:end))
kx=linspace(kpoints(1),kpoints(end),1000);
h=plot(kx,dfunc(an,kx),'r-','linewidth',2);
hold off
xlim([0.42,0.5])
%ylim([0 800])
set(gca,'fontsize',12)
legend('MPB group velocity data','Curve fit','location','best')
text(0.433,550,'$f=1/(2a_{1}(x-0.5)+4a_{2}(x-0.5)^{3}+6a_{3}(x-0.5)^{5})$','Interpreter','latex', ...
    'fontsize',11)
text(0.433,500,'$a_{1}=0.0394,a_{2}=-106.3288,a_{3}=6.5887e3$','Interpreter','latex','fontsize',12)
saveas(h,'GVfit.png')

%% Effective mass

figure(2)
h1=plot(kpoints,data.bands(:,20),'o')
hold on
h2=plot(kpoints,data.bands(:,21))
h3=plot(kpoints,data.bands(:,22))
h5=area(kpoints,data.bands(:,19))
xx=[kpoints',fliplr(kpoints')];
inBetween1=[data.vectors(:,1)',fliplr(repmat(0.6,size(kpoints,1),1)')];
inBetween2=[data.bands(:,23)',fliplr(repmat(0.6,size(kpoints,1),1)')];
h9=fill(xx,inBetween2,[0.8 0.8 0.8],'EdgeColor','none')
h10=fill(xx,inBetween1,[1 1 1],'EdgeColor','none')
h6=plot(kpoints,data.vectors(:,1))
xlim([0.42 0.5])
ylim([0.39 0.5])

func=@(c,x)c(1)+(an(1)).*(x-0.5).^2+(an(2)).*(x-0.5).^4+(an(3)).*(x-0.5).^6;
c0=[1];
bn=lsqcurvefit(func,c0,kpoints(end-15:end),data.bands(end-15:end,20))
kx=linspace(kpoints(1),kpoints(end),1000);
h=plot(kx,func(bn,kx),'r-','linewidth',2);

hold off

dk=kx(1000)-kx(999);
dw=diff(func(bn,kx));
d2w=diff(dw./dk);
alpha=d2w./dk;
alpha(end)
SIunit_wb_even=gbands20(end).*SIunit_omega;
inv_eff_mass=alpha.*SIunit_wb_even/(SIunit_k0.^2);
delta=(inv_eff_mass(end)/(a^2)); %Units in THz
inv_eff_mass(end);
D=delta(end)*10^(-12);

set(h1,'Color',[1 0.2 0.1],'Linewidth',2)
set(h2,'Color',[0.6 0.6 0.6],'Linewidth',1.5)
set(h3,'Color',[0.1 0.4 1],'Linewidth',2)
set(h6,'Color',[0 0 0],'Linewidth',1.5)
set(h5,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
set(gca,'Layer','top')
set(gca,'Fontsize',12)
hText = text(0.425, 0.4426,'$$\nu_{D2}$$','Interpreter','Latex','FontSize',16);
text(0.46,0.42,'$$ L=\sqrt{\alpha \omega_{b}/k_{0}^2\Delta_{e}} $$','Interpreter','latex','fontsize',16)
text(0.46,0.41,['$$\Delta_{e}=2\pi\times$$' num2str(D,2) 'THz'],'Interpreter','latex','fontsize',16)
xlabel('Propagation constant ($k_{x}a_{0}/2\pi$)','Interpreter','Latex','FontSize',16)
ylabel('Normalized Frequency ($\omega a{_0}/2\pi c$)','Interpreter','Latex','FontSize',16)

set(gcf, 'PaperSize', [4 2]);
posLabelx=get(get(gca,'xlabel'),'Position')
posLabely=get(get(gca,'ylabel'),'Position')
set(get(gca,'xlabel'),'Position',posLabelx-[0,0.001,0])
set(get(gca,'ylabel'),'Position',posLabely-[0.0005,0,0])
set(gca,'Position',[0.15 0.15 0.7 0.8])
saveas(h1,'dispersionfit.png')

%% Quality factor estimation
numcell=30;
k0=0.5*(numcell-1)/numcell;
[M index]=min(abs(kx-k0));
vg=dfunc(an,kx(index)).^(-1);
wg=func(bn,kx(index))*SIunit_omega;
ng=1/vg;
R=(ng-1)/(ng+1);
L=12*10^-6;
Qfactor=(-2*ng*L*wg/c)*(1/log(R^2));