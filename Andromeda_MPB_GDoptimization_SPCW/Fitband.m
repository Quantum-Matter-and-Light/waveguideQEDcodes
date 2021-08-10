function [Delta,alpha,Meff]=Fitband(nband,nkpoints,latticeC)
%% Constant and Unit conversion
c=299792458;
% latticeC=374*10^-9;     % Lattice constant in m
SIunit_omega=(c/latticeC)*2*pi;  %%% frequency in THz
SIunit_k0=0.5*2*pi/latticeC;

%% Load MAT file
% nband=20;
% nkpoints=16;

load('data.mat');
data=data{1};
lightcone=data.vectors(:,1);
kpoints=data.vectors(:,1);
%D2line=repmat(d2_freq,1,length(kpoints));
band=data.bands(:,nband);
G_index=(data.xvelocity(end-nkpoints:end-1,nband)).^(-1);
kx=linspace(kpoints(1),kpoints(end),1000);

%% Polynomial Fit function (Three steps.) 
%% 1.fit dispersion(get param1). 
Nint=10;
eqstr='';
for i1=1:Nint
    if i1<Nint
        eqstr=[eqstr 'a(' num2str(i1) ').*(x-0.5).^(' num2str(2*i1) ')+'];
    else
        eqstr=[eqstr 'a(' num2str(i1) ').*(x-0.5).^(' num2str(2*i1) ')'];
    end
end
eqstr=['@(a,x)' eqstr '+a(' num2str(Nint+1) ')'];
Bfun=str2func(eqstr);
a0=repmat(-10,Nint+1,1);
[param1, resnorm, ~, exitflag, lsqoutput]=lsqcurvefit(Bfun,a0,kpoints(end-nkpoints:end-1),band(end-nkpoints:end-1))
param1_results{i1}=param1;

A=figure(1)
h=plot(kpoints,band,'bo','linewidth',2);
hold on
h=plot(kx,Bfun(param1,kx),'r-','linewidth',2);
hold off
grid on
title('Dispersion Fitting for Initial parameter','Interpreter','Latex','Fontsize',12);
xlabel('Propagation constant ($k_{x}a_{0}/2\pi$)','Interpreter','Latex','FontSize',16);
ylabel(['Normalized Frequency $$(\omega a_{0}/2 \pi c)$$'],'Interpreter','Latex','FontSize',16);
legend('MPB Dispersion Data','Curve fit','Location','Best');
saveas(A,'InitialFit.png');

%% 2.use param1 for initial value to fit group index (get param2.)
deqstr='';
for i1=1:Nint
    if i1<Nint
        deqstr=[deqstr '(' num2str(2*i1) ').*a(' num2str(i1) ').*(x-0.5).^(' num2str(2*i1-1) ')+'];
    else
        deqstr=[deqstr '(' num2str(2*i1) ').*a(' num2str(i1) ').*(x-0.5).^(' num2str(2*i1-1) ')'];
    end
end
deqstr=['@(a,x)(' deqstr ').^(-1)'];
dfun=str2func(deqstr);
[param2, resnorm, ~, exitflag, lsqoutput]=lsqcurvefit(dfun,param1,kpoints(end-nkpoints:end-1),G_index)
param2_results{i1}=param2;
kx=linspace(kpoints(1),kpoints(end),1000);

B=figure(2)
h=plot(kpoints(end-nkpoints:end-1),G_index,'bo','linewidth',2)
xlabel('Propagation constant ($k_{x}a_{0}/2\pi$)','Interpreter','Latex','FontSize',16)
ylabel('$1/v_{g}$','Interpreter','Latex','FontSize',16)
hold on
h=plot(kx,dfun(param2,kx),'r-','linewidth',1.5);
hold off
set(gca,'fontsize',12)
legend('MPB Group Velocity Data','Curve fit','location','best')
grid on
saveas(B,'GVfit.png')

%% 3. fit dispersion again with param2.)
Neqstr='';
for i1=1:Nint
    if i1<Nint
        Neqstr=[Neqstr '(' num2str(param2(i1)) ').*(x-0.5).^(' num2str(2*i1) ')+'];
    else
        Neqstr=[Neqstr '(' num2str(param2(i1)) ').*(x-0.5).^(' num2str(2*i1) ')'];
    end
end

b0=-0.1;
Neqstr=['@(b,x)' Neqstr '+b(' num2str(Nint+1) ')'];
Nfun=str2func(Neqstr);
bn=lsqcurvefit(Nfun,a0,kpoints(end-nkpoints:end-1),band(end-nkpoints:end-1));

%% Effective mass evaluation
dk=kx(1000)-kx(999);
dw=diff(Nfun(bn,kx));
d2w=diff(dw./dk);
alpha=d2w./dk;
alpha=alpha(end);
SIunit_wb=band(end).*SIunit_omega;
inv_eff_mass=alpha.*SIunit_wb/(SIunit_k0.^2);
Meff=inv_eff_mass.^(-1);
delta=(inv_eff_mass(end)/(latticeC^2)); %Units in THz
inv_eff_mass(end);
Delta=delta(end)*10^(-12);
ylim=get(gca,'ylim');
xlim=get(gca,'xlim')
grid on
set(gca,'Layer','top')
set(gca,'Fontsize',12)

%% Plot
C=figure(3)
h=plot(kpoints,band,'ro','linewidth',2)
hold on
h=plot(kx,Nfun(bn,kx),'b-','linewidth',2)
hold off
grid on
title('Dispersion Final Fitting','Interpreter','Latex','Fontsize',12);
xlabel('Propagation constant ($k_{x}a_{0}/2\pi$)','Interpreter','Latex','FontSize',16);
ylabel(['Normalized Frequency $$(\omega a_{0}/2 \pi c)$$'],'Interpreter','Latex','FontSize',16);
legend('MPB Dispersion Data','Curve fit');
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1)+0.01,ylim(2)-0.001,'$$ L=\sqrt{\alpha \omega_{b}/k_{0}^2\Delta_{e}} $$','Interpreter','latex','fontsize',16)
text(xlim(1)+0.01,ylim(2)-0.0015,['$$\Delta_{e}=2\pi\times$$' num2str(Delta,3) 'THz'],'Interpreter','latex','fontsize',16)
saveas(C,'FinalFit.png');

%% Quality factor estimation
numcell=40;
k0=0.5*(numcell-1)/numcell;
[M index]=min(abs(kx-k0));
vg=dfun(param2,kx(index)).^(-1);
wg=Nfun(bn,kx(index))*SIunit_omega;
ng=1/vg;
R=(ng-1)/(ng+1);
L=12*10^-6;
Qfactor=(-2*ng*L*wg/c)*(1/log(R^2));
