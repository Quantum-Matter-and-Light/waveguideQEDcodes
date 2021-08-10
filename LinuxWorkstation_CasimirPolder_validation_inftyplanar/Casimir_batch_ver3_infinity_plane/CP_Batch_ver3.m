%% 2017.08.22 Youn Seok Lee : Casimir-Polder potential calculation.

%% Directory and CTL file %%
workname='infty_plane';
fname_xyr='infty_plane';
% fname_yzr='Squircle_W1_yzr';
currdir=pwd;
xyctlname=strcat(fname_xyr,'.ctl');
% yzctlname=strcat(fname_yzr,'.ctl');
ctldir='/home/uqml/Desktop/MEEP/infty_plane';    %% CTL file location
workspace=strcat(ctldir,'/20170926_testend3Courant0.5res32_',workname);
outDielecdir_xyr=strcat(workspace,'/out_xyr');
outVacdir_xyr=strcat(workspace,'/out_vac_xyr');
% outDielecdir_yzr=strcat(workspace,'/out_yzr');
% outVacdir_yzr=strcat(workspace,'/out_vac_yzr');

if ~isdir(workspace)
    mkdir(workspace)
end

unix(['cp ' currdir '/spdf.m ' workspace '/spdf.m'])
unix(['cp ' currdir '/pshalf.mat ' workspace '/pshalf.mat'])
unix(['cp ' currdir '/qoperators.m ' workspace '/qoperators.m'])
unix(['cp ' currdir '/lgwt.m ' workspace '/lgwt.m'])
unix(['cp ' currdir '/Wigner6jcoeff.m ' workspace '/Wigner6jcoeff.m'])
unix(['cp ' currdir '/CP_groundstate_shift.m ' workspace '/CP_groundstate_shift.m'])
unix(['cp ' currdir '/vectortensor_Casimir.m ' workspace '/vectortensor_Casimir.m'])
unix(['cp ' currdir '/readEfield.m ' workspace '/readEfield.m'])
unix(['cp ' currdir '/readEfieldNk.m ' workspace '/readEfieldNk.m'])
unix(['cp ' currdir '/makegt.m ' workspace '/makegt.m'])
cd(ctldir)
unix(['cp ' ctldir '/' fname_xyr '.ctl ' workspace '/' fname_xyr '.ctl'])
% unix(['cp ' ctldir '/' fname_yzr '.ctl ' workspace '/' fname_yzr '.ctl'])
cd(workspace)

%% Constants %%
ep0=8.85e-12;
mu0=4*pi*10^(-7);
hbar=1.05457148e-34;
c=299792458;                                %Speed of light
kB=1.38064852*10^(-23);                     %Boltzmann constant
h=6.626e-34;                                %Planck's constant
pshalf=importdata('pshalf.mat');            %wavelength and life time for D1 and D2 line
p_index=pshalf(:,1);                        %principle quantum number
omega_ps=(2*pi*c)./(pshalf(:,[2 4])*1e-9);  %transition frequencies for D1 (first column) and D2 (second column) line
tau_ps=pshalf(:,[3,5])*1e-6;                %lifetime for D1 and D2 line
gamma_ps=1./tau_ps;                         %spontaneous emission rates

%% Computation parameter setting %%
a_meep=373*10^-9;                                %Lattice constant
N=1;                                       %Number of k-points
ka=-0.5;                                    
kb=0.5;
[kx,w]=lgwt(N,ka,kb);                       %Gaussian-quadrature integration
endtime=10;
res=32;
Courant=0.25;
dt=Courant/res;
t=0:dt:endtime;
xbgn=0%:(1/res):0.19/(a_meep*10^6);
ybgn=0%:(1/res):0.1/(a_meep*10^6);
zbgn=0%;
Sigma=1;
% rminy=4/res;
% rmaxy=36/res;
% dry=4/res;
%position_x=0:(2/res):0.2/0.335; % scanning range in x-direction
position_y=4/res:(8/res):1/(a_meep*10^6); %                   y-direction
position_z=0:(1/res):0.4/(a_meep*10^6); %                   z-direction

%% Unix commend: running MEEP %%
np=6;                                       % number of cores allocated for each run

for ii=1:1%length(kx)/2
    for jj=1:length(xbgn)
    unix(['mpirun -np ' num2str(np) ' meep-openmpi kpt=' num2str(kx(ii)) ' res=' num2str(res) ...
        ' r-min=' num2str(position_y(1)) ' r-max=' num2str(position_y(end)) ' dr=' num2str(8/res) ' xbgn=' num2str(xbgn(jj)) ' endtime=' num2str(endtime) ...
        ' crnt=' num2str(Courant) ' Sigma=' num2str(Sigma) ' ' xyctlname ' |tee k' num2str(kx(ii)) '_xr' num2str(xbgn(jj)) '_' fname_xyr '.out'])
    end
%     for jj=1:length(ybgn)
%     unix(['mpirun -np ' num2str(np) ' meep-mpi kpt=' num2str(kx(ii)) ' res=' num2str(res) ...
%         ' r-min=' num2str(position_z(1)) ' r-max=' num2str(position_z(end)) ' dr=' num2str(1/res) ' ybgn=' num2str(ybgn(jj)) ' endtime=' num2str(endtime) ...
%         ' crnt=' num2str(Courant) ' Sigma=' num2str(Sigma) ' ' yzctlname ' |tee k' num2str(kx(ii)) '_yr' num2str(ybgn(jj)) '_' fname_yzr '.out'])
%     end
end 

%% Output file compilation // Data import // Gaussian quadrature integration %%

for ii=1:1%length(kx)/2
    for kk=1:length(xbgn)
        if kx==0  % non-periodic structure (e.g. infinite plane)
            unix(['cp ' workspace '/readEfieldNk.m ' strcat(outDielecdir_xyr,'_',num2str(xbgn(kk))) '/readEfieldNk.m'])
            unix(['cp ' workspace '/readEfieldNk.m ' strcat(outVacdir_xyr,'_',num2str(xbgn(kk))) '/readEfieldNk.m'])
            cd(strcat(outDielecdir_xyr,'_',num2str(xbgn(kk))))
            [xyDielecEtensor]=readEfieldNk(length(position_y),kx(ii));
            DielecEtensor_xyr(:,:,:,:,kk,ii)=xyDielecEtensor(:,:,:,:).*(mu0*c^2/a_meep^3);
            cd(strcat(outVacdir_xyr,'_',num2str(xbgn(kk))))
            [xyVacEtensor]=readEfieldNk(length(position_y),kx(ii));
            VacEtensor_xyr(:,:,:,:,kk,ii)=xyVacEtensor(:,:,:,:).*(mu0*c^2/a_meep^3);
        else
            unix(['cp ' workspace '/readEfield.m ' strcat(outDielecdir_xyr,'_',num2str(xbgn(kk))) '/readEfield.m'])
            unix(['cp ' workspace '/readEfield.m ' strcat(outVacdir_xyr,'_',num2str(xbgn(kk))) '/readEfield.m'])
            cd(strcat(outDielecdir_xyr,'_',num2str(xbgn(kk))))
            [xyDielecEtensor]=readEfield(length(position_y),kx(ii),w(ii));
            DielecEtensor_xy(:,:,:,:,kk,ii)=xyDielecEtensor(:,:,:,:);
            cd(strcat(outVacdir_xyr,'_',num2str(xbgn(kk))))
            [xyVacEtensor]=readEfield(length(position_y),kx(ii),w(ii));
            VacEtensor_xy(:,:,:,:,kk,ii)=xyVacEtensor(:,:,:,:);
            
        end
    end
%     for jj=1:length(ybgn)
%         if kx==0  % non-periodic structure (e.g. infinite plane)
%             unix(['cp ' workspace '/readEfieldNk.m ' strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))) '/readEfieldNk.m'])
%             unix(['cp ' workspace '/readEfieldNk.m ' strcat(outVacdir_yzr,'_',num2str(ybgn(jj))) '/readEfieldNk.m'])
%             cd(strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))))   
%             [yzDielecEtensor]=readEfieldNk(length(position_z),kx(ii));
%             DielecEtensor_yzr(:,:,:,:,jj,ii)=yzDielecEtensor(:,:,:,:);
%             cd(strcat(outVacdir_yzr,'_',num2str(ybgn(jj))))    
%             [yzVacEtensor]=readEfieldNk(length(position_z),kx(ii));
%             VacEtensor_yzr(:,:,:,:,jj,ii)=yzVacEtensor(:,:,:,:);
%         else
%             unix(['cp ' workspace '/readEfield.m ' strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))) '/readEfield.m'])
%             unix(['cp ' workspace '/readEfield.m ' strcat(outVacdir_yzr,'_',num2str(ybgn(jj))) '/readEfield.m'])
%             cd(strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))))   
%             [yzDielecEtensor]=readEfield(length(position_z),kx(ii),w(ii));
%             DielecEtensor_yz(:,:,:,:,jj,ii)=yzDielecEtensor(:,:,:,:);
%             cd(strcat(outVacdir_yzr,'_',num2str(ybgn(jj))))    
%             [yzVacEtensor]=readEfield(length(position_z),kx(ii),w(ii));
%             VacEtensor_yz(:,:,:,:,jj,ii)=yzVacEtensor(:,:,:,:);
%         end
%     end
end

%clear xyDielecEtensor xyVacEtensor yzDielecEtensor yzVacEtensor
if kx~=0
    for i=1:size(DielecEtensor_yz,1)
        for j=1:size(DielecEtensor_yz,2)
            for k=1:size(DielecEtensor_yz,3)
                for l=1:size(DielecEtensor_xy,4)
                    for m=1:size(DielecEtensor_xy,5)
                        DielecEtensor_xyr(i,j,k,l,m)=2*sum(DielecEtensor_xy(i,j,k,l,m,:).*(mu0*c^2/a_meep^3));
                        VacEtensor_xyr(i,j,k,l,m)=2*sum(VacEtensor_xy(i,j,k,l,m,:).*(mu0*c^2/a_meep^3));
                    end
                end;
%                 for l=1:size(DielecEtensor_yz,4)
%                     for m=1:size(DielecEtensor_yz,5)
%                         DielecEtensor_yzr(i,j,k,l,m)=2*sum(DielecEtensor_yz(i,j,k,l,m,:).*(mu0*c^2/a_meep^3));
%                         VacEtensor_yzr(i,j,k,l,m)=2*sum(VacEtensor_yz(i,j,k,l,m,:).*(mu0*c^2/a_meep^3));
%                     end
%                 end;
            end;
        end;
    end;
end
%% Dipole reduced matrix elements for each transitions %%

cd(workspace)
[Fi,mfi,omega_ff,Fi_easy,p_index,FF,alpha1_coeff,alpha2_coeff]=CP_groundstate_shift();
alpha1_coeff=alpha1_coeff*(-1);
for ii=1:6                                  %index(ii): principle quantum number 
    for ind_ff=1:2                              %index(jj): ground state hyperfine splitting (F=3 or 4)
        alpha0_1(:,:,ii)=(abs(FF(:,:,ii)).^2)*(2/(3));
        alpha0(ii,ind_ff)=sum(alpha0_1(ind_ff,1:2,ii));
        alpha0(ii+6,ind_ff)=sum(alpha0_1(ind_ff,3:6,ii));
           
        alpha1_1(:,:,ii)=((alpha1_coeff(:,:,ii).*abs(FF(:,:,ii)).^2));
        alpha1(ii,ind_ff)=sum(alpha1_1(ind_ff,1:2,ii));
        alpha1(ii+6,ind_ff)=sum(alpha1_1(ind_ff,3:6,ii));
    
        alpha2_1(:,:,ii)=((alpha2_coeff(:,:,ii).*abs(FF(:,:,ii)).^2));
        alpha2(ii,ind_ff)=sum(alpha2_1(ind_ff,1:2,ii));
        alpha2(ii+6,ind_ff)=sum(alpha2_1(ind_ff,3:6,ii));
    end
end

%% Construction of g(-t) function %%
omega_ps_meep=omega_ps.*a_meep./c;          %Resonant frequncies in MEEP unit
omega_ps_meep(7:12,1)=omega_ps_meep(:,2);
omega_trans(:,1)=omega_ps_meep(:,1);
for kk=1:12
    trans_freq=omega_trans(kk);
    [g]=makegt(Courant,res,endtime,trans_freq,Sigma,a_meep);
    
    %% Gamma_ij construction %%
    for i=1:3
        for j=1:3
            for l=1:length(position_y)
                for m=1:length(xbgn)
                    Gamma_xy(i,j,l,m,kk)=sum(imag(g(:)).*real(DielecEtensor_xyr(:,i,j,l,m)).*dt.*a_meep/c);
                    Gamma_vac_xy(i,j,l,m,kk)=sum(imag(g(:)).*real(VacEtensor_xyr(:,i,j,l,m)).*dt.*a_meep/c);
                end
            end
%             for p=1:length(position_z)
%                 for q=1:length(ybgn)
%                     Gamma_yz(i,j,p,q,kk)=sum(imag(g(:)).*real(DielecEtensor_yzr(:,i,j,p,q)).*dt.*a_meep/c);
%                     Gamma_vac_yz(i,j,p,q,kk)=sum(imag(g(:)).*real(VacEtensor_yzr(:,i,j,p,q)).*dt.*a_meep/c);
%                 end
%             end
        end
    end
    
    %% CP potential (Spherical tensor decomposition)%%
    for ind_ff=1:2
        mfI=mfi{1};
        mf_i1=mfI{ind_ff}';
        for i=1:length(position_y)
            for j=1:length(xbgn)
                Green_xyr=Gamma_xy(:,:,i,j,kk);
                Vacuum_xyr=Gamma_vac_xy(:,:,i,j,kk);
                % Green_xyr=Gamma_xy(:,:,i,kk);
                % Vacuum_xyr=Gamma_vac_xy(:,:,i,kk);
                if ind_ff==1
                states1=[repmat(Fi(ind_ff),numel(mf_i1),1),mf_i1];
                [scalar_coeff_1,vector_coeff_1,tensor_coeff_1]=vectortensor_Casimir(states1,Fi(1),Green_xyr);
                [scalar_coeff_2,vector_coeff_2,tensor_coeff_2]=vectortensor_Casimir(states1,Fi(1),Vacuum_xyr);
                E_shift_F_3_CP_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                E_shift_F_3_Vac_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
                elseif ind_ff==2
                states2=[repmat(Fi(ind_ff),numel(mf_i1),1),mf_i1];
                [scalar_coeff_1,vector_coeff_1,tensor_coeff_1]=vectortensor_Casimir(states2,Fi(2),Green_xyr);
                [scalar_coeff_2,vector_coeff_2,tensor_coeff_2]=vectortensor_Casimir(states2,Fi(2),Vacuum_xyr);
                E_shift_F_4_CP_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                E_shift_F_4_Vac_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
                end
            end
        end
%         for l=1:length(position_z)
%             for m=1:length(ybgn)
%                 Green_yzr=Gamma_yz(:,:,l,m,kk);
%                 Vacuum_yzr=Gamma_vac_yz(:,:,l,m,kk);
%                 if ind_ff==1
%                     states1=[repmat(Fi(ind_ff),numel(mf_i1),1),mf_i1];
%                     [scalar_coeff_1,vector_coeff_1,tensor_coeff_1]=vectortensor_Casimir(states1,Fi(1),Green_yzr);
%                     [scalar_coeff_2,vector_coeff_2,tensor_coeff_2]=vectortensor_Casimir(states1,Fi(1),Vacuum_yzr);
%                     E_shift_F_3_CP_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
%                     E_shift_F_3_Vac_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
%                 elseif ind_ff==2
%                     states2=[repmat(Fi(ind_ff),numel(mf_i1),1),mf_i1];
%                     [scalar_coeff_1,vector_coeff_1,tensor_coeff_1]=vectortensor_Casimir(states2,Fi(2),Green_yzr);
%                     [scalar_coeff_2,vector_coeff_2,tensor_coeff_2]=vectortensor_Casimir(states2,Fi(2),Vacuum_yzr);
%                     E_shift_F_4_CP_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
%                     E_shift_F_4_Vac_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
%                 end
%             end
%         end
    end
end

%clear DielecEtensor_xyr DielecEtensor_zr VacEtensor_xyr VacEtensor_zr
%clear Gamma_xy Gamma_vac_xy Gamma_z Gamma_vac_z Green_xyr Green_zr Vacuum_xyr Vacuum_zr

%% Optical Stark shift from trap mode by MPB field profile
% cd('/home/uqml/Desktop/MEEP/CP_final/Final_figure')
% amplitude=10*10^6;
% nband=22;
% kpoint=2;
% [HHXY_stark,HHYZ_stark]=Starkshift(xbgn,ybgn,position_y,position_z,amplitude,nband,kpoint,a_meep)
% H_stark_xy_F3=HHXY_stark{1,1};
% H_stark_xy_F4=HHXY_stark{1,2};
% H_stark_yz_F3=HHYZ_stark{1,1};
% H_stark_yz_F4=HHYZ_stark{1,2};
% cd(workspace)
%% sum over all atomic transition // Tensor diagonlization %%

for l=1:length(position_y)
    for k=1:length(xbgn)
       for i=1:7
            for j=1:7
                E_total_shift_F_3_CP_xy(i,j,k,l)=sum(E_shift_F_3_CP_xy(i,j,l,k,:));
                E_total_shift_F_3_Vac_xy(i,j,k,l)=sum(E_shift_F_3_Vac_xy(i,j,l,k,:));
            end
        end
        for i=1:9
            for j=1:9
                E_total_shift_F_4_CP_xy(i,j,k,l)=sum(E_shift_F_4_CP_xy(i,j,l,k,:));
                E_total_shift_F_4_Vac_xy(i,j,k,l)=sum(E_shift_F_4_Vac_xy(i,j,l,k,:));
            end
        end
    [QCxy3,CP_xy_F3(:,:,k,l)]=eig(E_total_shift_F_3_CP_xy(:,:,k,l));
    [QVxy3,Vac_xy_F3(:,:,k,l)]=eig(E_total_shift_F_3_Vac_xy(:,:,k,l));
    [QCxy4,CP_xy_F4(:,:,k,l)]=eig(E_total_shift_F_4_CP_xy(:,:,k,l));
    [QVxy4,Vac_xy_F4(:,:,k,l)]=eig(E_total_shift_F_4_Vac_xy(:,:,k,l));
    CP_potential_xy(k,l,1)=trace(CP_xy_F3(:,:,k,l))-trace(Vac_xy_F3(:,:,k,l));
    CP_potential_xy(k,l,2)=trace(CP_xy_F4(:,:,k,l))-trace(Vac_xy_F4(:,:,k,l));
    Casimir_potential_xy=CP_potential_xy(:,:,1)+CP_potential_xy(:,:,2);
    
    end
end
save('CP_potential_xy.mat','Casimir_potential_xy')

%% Unit conversion // PLOT %%
CP_freq_unit_xy=real(Casimir_potential_xy)./h.*10.^(-3)/(2*pi); %units in KHz

xx=position_y.*(a_meep*10.^6);
figure(3)
plot(xx,CP_freq_unit_xy,'bo')
save('Sigma1.mat','CP_freq_unit_xy');
hold on
fun=@(c,x)(-1*c(2))./((x-c(3)).*x.^(3));
c0=[0 0 0];
c=lsqcurvefit(fun,c0,xx,CP_freq_unit_xy);
x=linspace(xx(1),xx(end),1000);
plot(x,fun(c,x),'r-','linewidth',1.5)
hold off
legend('Numercal data','curve fitting')
xlabel('d (um)')
ylabel('CP potential (kHz)')

f=c(3) 
C3=c(2)/c(3)    %1.16 kHz*um^3
C4=c(2) %0.15 kHz*um^4

figure(4)
x1=linspace(0.001,100,100000);
CP_planar=real(CP_freq_unit_xy).*xx.^3*-1;
func=@(C,X)(-1*C(2))./(X-C(3));
C0=[0 0 0];
C=lsqcurvefit(func,C0,xx,CP_planar);
loglog(xx,CP_planar,'ro')
hold on
loglog(x1,-1*C(2)./(x1-C(3)),'k-')
%hold off
legend('Numercal data','curve fitting')
xlabel('d (um)')
ylabel('Ucp*d^3 (kHz*um^3)')

CC3=C(2)/C(3)
CC4=C(2)
