%function [Gamma,Jij]=GammaJij_ver4(GN)
%% Constants
clight=299792458;
mu0=4*pi*10^(-7);
ep0=8.85e-12;
hbar=1.05457148e-34;
h=2*pi*hbar;
fcs_d1=335.116*10^12;
fcs_d2=351.725*10^12; % Units in Hz
omega_d1=2*pi*fcs_d1;
omega_d2=2*pi*fcs_d2;
ec=1.6022*10^-19;
a0=5.291*10^-11;

%% Oscillator strength for D1 & D2 transition including Zeeman sublevels.
I=7/2;
Li=0; %S
[Li,Ji,Fi,mfi]=spdf(Li);
Lf=1; %P
[Lf,Jf,Ff,mff]=spdf(Lf);
Fi=Fi{1}; % hyperfine ground states
Ffd1=Ff{1}; % hyperfine excited states for D1 transition
Ffd2=Ff{2}; % hyperfine excited states for D2 transition
mfig=mfi{1};
mff_d1=mff{1};
mff_d2=mff{2};

%JdJD1=4.489*ec*a0; % Jonathan's thesis (A.19)
%JdJD2=6.324*ec*a0; % Jonathan's thesis (A.20)

Ff_eas=[Ff{1},Ff{2}];
Fi_easy=repmat(Fi',1,numel(Ff_eas));
Ff_easy=repmat(Ff_eas,numel(Fi),1);
Ji_easy=repmat([Ji;Ji],1,numel(Ff_eas));
Jf_easy=repmat([Jf(1)*ones(size(Ff{1})),Jf(2)*ones(size(Ff{2}))],numel(Jf),1);

pshalf=importdata('pshalf.mat');
omega_ps=(2*pi*clight)./(pshalf(1,[2 4])*1e-9);
tau_ps=pshalf(1,[3,5])*1e-6;
gamma_ps=1./tau_ps;
jj_factor=(2*Ji+1)./(2*Jf(2)+1);
jj_ps_sq=(3*pi*ep0*hbar*clight^3)./(omega_ps.^3).*(1./(jj_factor)).*gamma_ps;
Dm_d1=sqrt(jj_ps_sq(1,1));
Dm_d2=sqrt(jj_ps_sq(1,2));
q=[-1 0 1];

%% Oscillator strength
for ind_i=1:numel(Fi)
    for ind_f=1:numel(Ffd1)
        mfi1=mfig{ind_i};
        mffd1=mff_d1{ind_f};
        Hd1=zeros(numel(mfi1),numel(mffd1));
        HdH1=zeros(numel(mfi1),numel(mffd1));
        for ii=1:3
            for mg=1:numel(mfi1)
                for me=1:numel(mffd1)
                    [FdF(mg,me,ii)]=DipoleZeeman(0.5,0.5,Fi(ind_i),Ffd1(ind_f),mffd1(me),mfi1(mg),I,q(ii));
                    Hd1(mg,me,ii)=FdF(mg,me,ii);
                    HdH1(mg,me)=HdH1(mg,me)+FdF(mg,me,ii);
                    
                end
            end
        end
        dipole_d1_Zeeman{ind_i,ind_f}=Hd1;
        dipole_d1_hyperfine(ind_i,ind_f)=trace(HdH1'*HdH1)*(Dm_d1)^2;
    end
end

for ind_i=1:numel(Fi)
    for ind_f=1:numel(Ffd2)
        mfi1=mfig{ind_i};
        mffd2=mff_d2{ind_f};
        Hd2=zeros(numel(mfi1),numel(mffd2));
        HdH2=zeros(numel(mfi1),numel(mffd2));
        for ii=1:3
            for mg=1:numel(mfi1)
                for me=1:numel(mffd2)
                    [FdF(mg,me,ii)]=DipoleZeeman(0.5,1.5,Fi(ind_i),Ffd2(ind_f),mffd2(me),mfi1(mg),I,q(ii));
                    Hd2(mg,me,ii)=FdF(mg,me,ii);
                    HdH2(mg,me)=HdH2(mg,me)+FdF(mg,me,ii);
                    
                end
            end
        end
        dipole_d2_Zeeman{ind_i,ind_f}=Hd2;
        dipole_d2_hyperfine(ind_i,ind_f)=trace(HdH2'*HdH2)*(Dm_d2)^2;        
    end
end




%% Total Decay rate and Lamb shift for 5P3/2 state
% G=omega_d2/(6*pi*clight)*eye(3);
%  Gamma=sum(sum(dipole_d2_hyperfine))*(2*omega_d2^2/(hbar*ep0*clight^2))*(imag(G))*jj_factor;
%  Jij=sum(sum(dipole_d2_hyperfine))*(omega_d2^2/(hbar*ep0*clight^2))*(real(G))*jj_factor;

%% Decay rate and Lamb shift via Spherical decomposition and sum over all q,q'
GN=omega_d2/(6*pi*clight)*eye(3);
G0=omega_d2/(6*pi*clight)*eye(3);

GN00=(-1/sqrt(3))*trace(GN);
GN1p1=(1/2)*(GN(3,1)-GN(1,3)+1i*(GN(3,2)-GN(2,3)));
GN1m1=(1/2)*(GN(3,1)-GN(1,3)-1i*(GN(3,2)-GN(2,3)));
GN1z0=(1i/sqrt(2))*(GN(1,2)-GN(2,1));
GN2p2=(1/2)*(GN(1,1)-GN(2,2)+1i*(GN(1,2)+GN(2,1)));
GN2m2=(1/2)*(GN(1,1)-GN(2,2)-1i*(GN(1,2)+GN(2,1)));
GN2p1=(-1/2)*(GN(1,3)+GN(3,1)+1i*(GN(2,3)+GN(3,2)));
GN2m1=(1/2)*(GN(1,3)+GN(3,1)-1i*(GN(2,3)+GN(3,2)));
GN2z0=(1/sqrt(6))*(3*GN(3,3)-trace(GN));

for ind_i=1:numel(Fi)
    for ind_f=1:numel(Ffd2)
        temp=dipole_d2_Zeeman{ind_i,ind_f};
        dipGN_00=(-1/sqrt(3)).*(conj(temp(:,:,2))'*temp(:,:,2)+conj(temp(:,:,1))'*temp(:,:,1)+conj(temp(:,:,3))'*temp(:,:,3)).*GN00;
        dipGN_1z0=(1/sqrt(2)).*(conj(temp(:,:,3))'*temp(:,:,3)-conj(temp(:,:,1))'*temp(:,:,1)).*GN1z0.*(-1);
        dipGN_1p1=(-1/sqrt(2)).*(temp(:,:,2)'*conj(temp(:,:,3))+conj(temp(:,:,2))'*temp(:,:,3)).*GN1m1.*(-1);
        dipGN_1m1=(1/sqrt(2)).*(temp(:,:,2)'*conj(temp(:,:,1))+conj(temp(:,:,2))'*temp(:,:,1)).*GN1p1.*(-1);
        dipGN_2z0=(1/sqrt(6)).*(2*conj(temp(:,:,2))'*temp(:,:,2)-conj(temp(:,:,3))'*temp(:,:,3)-conj(temp(:,:,1))'*temp(:,:,1)).*GN2z0;
        dipGN_2p1=(-1/sqrt(2)).*(temp(:,:,2)'*conj(temp(:,:,3))-conj(temp(:,:,2))'*temp(:,:,3)).*GN2m1;
        dipGN_2m1=(-1/sqrt(2)).*(temp(:,:,2)'*conj(temp(:,:,1))-conj(temp(:,:,2))'*temp(:,:,1)).*GN2p1;
        dipGN_2p2=-1*temp(:,:,3)'*conj(temp(:,:,3)).*GN2m2;
        dipGN_2m2=-1*temp(:,:,1)'*conj(temp(:,:,1)).*GN2p2;
        dipGN=dipGN_00+dipGN_1z0+dipGN_1p1+dipGN_1m1+dipGN_2z0+dipGN_2p1+dipGN_2m1+dipGN_2p2+dipGN_2m2;
        GammaJijN{ind_i,ind_f}=dipGN.*(2*omega_d2^2/(hbar*ep0*clight^2))*(Dm_d2)^2*jj_factor;
        [Q,H]=eig(GammaJijN{ind_i,ind_f});
        GammaFF{ind_i,ind_f}=sort(imag(diag(H)));
        JijFF{ind_i,ind_f}=sort(real(diag(H)));
        GammaN(ind_i,ind_f)=sum(imag(diag(H)));
        JijN(ind_i,ind_f)=sum(real(diag(H)))/2;
    end
end
Gamma=sum(sum(GammaN))
Jij=sum(sum(JijN))
% Gamma=GammaFF{2,4}
% Jij=JijFF{2,4}
%% Decay rate and Lamb shift in Cartesian coordinates
% GN=omega_d2/(6*pi*clight)*eye(3);
% G0=omega_d2/(6*pi*clight)*eye(3);
% 
% for ind_i=1:numel(Fi)
%     for ind_f=1:numel(Ffd2)
%         temp=dipole_d2_Zeeman{ind_i,ind_f};
%         D2dx=(1/sqrt(2))*(temp(:,:,1)-temp(:,:,3));
%         D2dy=(1i/sqrt(2))*(temp(:,:,1)+temp(:,:,3));
%         D2dz=temp(:,:,2);
%         D2d=temp(:,:,1)+temp(:,:,2)+temp(:,:,3);
%         D2ddx=conj(D2dx)';
%         D2ddy=conj(-D2dy)';
%         D2ddz=conj(D2dz)';
%         D2Dxx{ind_i,ind_f}=D2ddx*D2dx;
%         D2Dyy{ind_i,ind_f}=D2ddy*D2dy;
%         D2Dzz{ind_i,ind_f}=D2ddz*D2dz;
%         D2D{ind_i,ind_f}=D2d'*D2d;
%         ex=[1 0 0];
%         ey=[0 1 0];
%         ez=[0 0 1];
%         GammaJijNxx{ind_i,ind_f}=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*diag(D2Dxx{ind_i,ind_f}).*(ex*GN*ex');
%         GammaJij0xx{ind_i,ind_f}=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*diag(D2Dxx{ind_i,ind_f}).*(ex*G0*ex');
%         GammaJijNyy{ind_i,ind_f}=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*diag(D2Dyy{ind_i,ind_f}).*(ey*GN*ey');
%         GammaJij0yy{ind_i,ind_f}=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*diag(D2Dyy{ind_i,ind_f}).*(ey*G0*ey');
%         GammaJijNzz{ind_i,ind_f}=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*diag(D2Dzz{ind_i,ind_f}).*(ez*GN*ez');
%         GammaJij0zz{ind_i,ind_f}=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*diag(D2Dzz{ind_i,ind_f}).*(ez*G0*ez');
%         GammaJijNHxx(ind_i,ind_f)=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*trace(D2D{ind_i,ind_f}).*(ex*GN*ex');
%         GammaJij0Hxx(ind_i,ind_f)=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*trace(D2D{ind_i,ind_f}).*(ex*G0*ex');
%         GammaJijNHyy(ind_i,ind_f)=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*trace(D2D{ind_i,ind_f}).*(ey*GN*ey');
%         GammaJij0Hyy(ind_i,ind_f)=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*trace(D2D{ind_i,ind_f}).*(ey*G0*ey');
%         GammaJijNHzz(ind_i,ind_f)=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*trace(D2D{ind_i,ind_f}).*(ez*GN*ez');
%         GammaJij0Hzz(ind_i,ind_f)=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*trace(D2D{ind_i,ind_f}).*(ez*G0*ez');
%     end
% end
% TGJxx=sum(GammaJij0xx{1,1})+sum(GammaJij0xx{1,2})+sum(GammaJij0xx{1,3})+sum(GammaJij0xx{1,4})+ ...
%      sum(GammaJij0xx{2,1})+sum(GammaJij0xx{2,2})+sum(GammaJij0xx{2,3})+sum(GammaJij0xx{2,4});
% TGJyy=sum(GammaJij0xx{1,1})+sum(GammaJij0yy{1,2})+sum(GammaJij0yy{1,3})+sum(GammaJij0yy{1,4})+ ...
%      sum(GammaJij0xx{2,1})+sum(GammaJij0yy{2,2})+sum(GammaJij0yy{2,3})+sum(GammaJij0yy{2,4});
% TGJzz=sum(GammaJij0zz{1,1})+sum(GammaJij0zz{1,2})+sum(GammaJij0zz{1,3})+sum(GammaJij0zz{1,4})+ ...
%      sum(GammaJij0zz{2,1})+sum(GammaJij0zz{2,2})+sum(GammaJij0zz{2,3})+sum(GammaJij0zz{2,4});

% eq=[(-1/sqrt(2))*(ex+1i*ey), (1/sqrt(2))*(ex-1i*ey), ez];
% GammaN=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*D2*imag(GN);
% Gamma0=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*D2*imag(G0);
%JijN=(omega_d2^2/(hbar*ep0*clight^2))*jj_factor*D2*real(GN);
%Jij0=(omega_d2^2/(hbar*ep0*clight^2))*jj_factor*D2*real(G0);
%GammaJijN=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*D2*(GN);
%GammaJij0=(2*omega_d2^2/(hbar*ep0*clight^2))*jj_factor*D2*(G0);

%% Cartesian and Spherical tensor expression
% GN00=(-1/sqrt(3))*trace(GN);
% GN1p1=(1/2)*(GN(3,1)-GN(1,3)+1i*(GN(3,2)-GN(2,3)));
% GN1m1=(1/2)*(GN(3,1)-GN(1,3)-1i*(GN(3,2)-GN(2,3)));
% GN1z0=(1i/sqrt(2))*(GN(1,2)-GN(2,1));
% GN2p2=(1/2)*(GN(1,1)-GN(2,2)+1i*(GN(1,2)+GN(2,1)));
% GN2m2=(1/2)*(GN(1,1)-GN(2,2)-1i*(GN(1,2)+GN(2,1)));
% GN2p1=(-1/2)*(GN(1,3)+GN(3,1)+1i*(GN(2,3)+GN(3,2)));
% GN2m1=(1/2)*(GN(1,3)+GN(3,1)-1i*(GN(2,3)+GN(3,2)));
% GN2z0=(1/sqrt(6))*(3*GN(3,3)-trace(GN));
% 
% G000=(-1/sqrt(3))*trace(G0);
% G01p1=(1/2)*(G0(3,1)-G0(1,3)+1i*(G0(3,2)-G0(2,3)));
% G01m1=(1/2)*(G0(3,1)-G0(1,3)-1i*(G0(3,2)-G0(2,3)));
% G01z0=(1i/sqrt(2))*(G0(1,2)-G0(2,1));
% G02p2=(1/2)*(G0(1,1)-G0(2,2)+1i*(G0(1,2)+G0(2,1)));
% G02m2=(1/2)*(G0(1,1)-G0(2,2)-1i*(G0(1,2)+G0(2,1)));
% G02p1=(-1/2)*(G0(1,3)+G0(3,1)+1i*(G0(2,3)+G0(3,2)));
% G02m1=(1/2)*(G0(1,3)+G0(3,1)-1i*(G0(2,3)+G0(3,2)));
% G02z0=(1/sqrt(6))*(3*G0(3,3)-trace(G0));
% 
% GNC(1,1)=(1/6)*((-2/sqrt(3))*GN00-sqrt(6)*GN2z0+3*(GN2m2+GN2p2));
% GNC(1,2)=(-1i/2)*(2*GN1z0-GN2m2+GN2p2);
% GNC(1,3)=(1/2)*(conj(GN1m1+GN1p1+GN2m1)-conj(GN2p1));
% GNC(2,1)=conj(GNC(1,2));
% GNC(2,2)=(1/6)*((-2/sqrt(3))*GN00-sqrt(6)*GN2z0-3*(GN2m2+GN2p2));
% GNC(2,3)=(-1i/2)*(-conj(GN1p1)+conj(GN1m1+GN2m1+GN2p1));
% GNC(3,1)=conj(GNC(1,3));
% GNC(3,2)=conj(GNC(2,3));
% GNC(3,3)=(1/sqrt(3))*(GN00+sqrt(2)*GN2z0);
% 
% G0C(1,1)=(1/6)*((-2/sqrt(3))*G000-sqrt(6)*GN2z0+3*(G02m2+G02p2));
% G0C(1,2)=(-1i/2)*(2*G01z0-G02m2+G02p2);
% G0C(1,3)=(1/2)*(conj(G01m1+G01p1+G02m1)-conj(G02p1));
% G0C(2,1)=conj(G0C(1,2));
% G0C(2,2)=(1/6)*((-2/sqrt(3))*G000-sqrt(6)*G02z0-3*(G02m2+G02p2));
% G0C(2,3)=(-1i/2)*(-conj(G01p1)+conj(G01m1+G02m1+G02p1));
% G0C(3,1)=conj(G0C(1,3));
% G0C(3,2)=conj(G0C(2,3));
% G0C(3,3)=(1/sqrt(3))*(G000+sqrt(2)*G02z0);


