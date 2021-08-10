%% 2017.08.22 Youn Seok Lee : Casimir-Polder potential calculation.

%% Directory and CTL file %%
workname='Squircle_W1';
fname_xyr='Squircle_W1_xyr';
fname_yzr='Squircle_W1_yzr';
currdir=pwd;
xyctlname=strcat(fname_xyr,'.ctl');
yzctlname=strcat(fname_yzr,'.ctl');
ctldir='/home/uqml/Desktop/MEEP/20170926_CP_final';    %% CTL file location
workspace=strcat(ctldir,'/20170917_end3Courant0.5res32_Finalrun_',workname);
outDielecdir_xyr=strcat(workspace,'/out_xyr');
outVacdir_xyr=strcat(workspace,'/out_vac_xyr');
outDielecdir_yzr=strcat(workspace,'/out_yzr');
outVacdir_yzr=strcat(workspace,'/out_vac_yzr');

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
unix(['cp ' ctldir '/' fname_yzr '.ctl ' workspace '/' fname_yzr '.ctl'])
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
N=10;                                       %Number of k-points
ka=-0.5;                                    
kb=0.5;
[kx,w]=lgwt(N,ka,kb);                       %Gaussian-quadrature integration
endtime=3;
res=32;
Courant=0.25;
dt=Courant/res;
t=0:dt:endtime;
xbgn=0:(1/res):0.19/(a_meep*10^6);
ybgn=0:(1/res):0.1/(a_meep*10^6);
zbgn=0;
Sigma=1;
% rminy=4/res;
% rmaxy=36/res;
% dry=4/res;
%position_x=0:(2/res):0.2/0.335; % scanning range in x-direction
position_y=0:(1/res):0.1/(a_meep*10^6); %                   y-direction
position_z=0:(1/res):0.4/(a_meep*10^6); %                   z-direction

%% Unix commend: running MEEP %%
np=16;                                       % number of cores allocated for each run

% for ii=1:length(kx)/2
%     for jj=1:length(xbgn)
%     unix(['mpirun -np ' num2str(np) ' meep-mpi kpt=' num2str(kx(ii)) ' res=' num2str(res) ...
%         ' r-min=' num2str(position_y(1)) ' r-max=' num2str(position_y(end)) ' dr=' num2str(1/res) ' xbgn=' num2str(xbgn(jj)) ' endtime=' num2str(endtime) ...
%         ' crnt=' num2str(Courant) ' Sigma=' num2str(Sigma) ' ' xyctlname ' |tee k' num2str(kx(ii)) '_xr' num2str(xbgn(jj)) '_' fname_xyr '.out'])
%     end
%     for jj=1:length(ybgn)
%     unix(['mpirun -np ' num2str(np) ' meep-mpi kpt=' num2str(kx(ii)) ' res=' num2str(res) ...
%         ' r-min=' num2str(position_z(1)) ' r-max=' num2str(position_z(end)) ' dr=' num2str(1/res) ' ybgn=' num2str(ybgn(jj)) ' endtime=' num2str(endtime) ...
%         ' crnt=' num2str(Courant) ' Sigma=' num2str(Sigma) ' ' yzctlname ' |tee k' num2str(kx(ii)) '_yr' num2str(ybgn(jj)) '_' fname_yzr '.out'])
%     end
% end 

%% Output file compilation // Data import // Gaussian quadrature integration %%

for ii=1:length(kx)/2
    for kk=1:length(xbgn)
        if kx==0  % non-periodic structure (e.g. infinite plane)
            unix(['cp ' workspace '/readEfieldNk.m ' strcat(outDielecdir_xyr,'_',num2str(xbgn(kk))) '/readEfieldNk.m'])
            unix(['cp ' workspace '/readEfieldNk.m ' strcat(outVacdir_xyr,'_',num2str(xbgn(kk))) '/readEfieldNk.m'])
            cd(strcat(outDielecdir_xyr,'_',num2str(xbgn(kk))))
            [xyDielecEtensor]=readEfieldNk(length(position_y),kx(ii));
            DielecEtensor_xyr(:,:,:,:,kk,ii)=xyDielecEtensor(:,:,:,:);
            cd(strcat(outVacdir_xyr,'_',num2str(xbgn(kk))))
            [xyVacEtensor]=readEfieldNk(length(position_y),kx(ii));
            VacEtensor_xyr(:,:,:,:,kk,ii)=xyVacEtensor(:,:,:,:);
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
    for jj=1:length(ybgn)
        if kx==0  % non-periodic structure (e.g. infinite plane)
            unix(['cp ' workspace '/readEfieldNk.m ' strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))) '/readEfieldNk.m'])
            unix(['cp ' workspace '/readEfieldNk.m ' strcat(outVacdir_yzr,'_',num2str(ybgn(jj))) '/readEfieldNk.m'])
            cd(strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))))   
            [yzDielecEtensor]=readEfieldNk(length(position_z),kx(ii));
            DielecEtensor_yzr(:,:,:,:,jj,ii)=yzDielecEtensor(:,:,:,:);
            cd(strcat(outVacdir_yzr,'_',num2str(ybgn(jj))))    
            [yzVacEtensor]=readEfieldNk(length(position_z),kx(ii));
            VacEtensor_yzr(:,:,:,:,jj,ii)=yzVacEtensor(:,:,:,:);
        else
            unix(['cp ' workspace '/readEfield.m ' strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))) '/readEfield.m'])
            unix(['cp ' workspace '/readEfield.m ' strcat(outVacdir_yzr,'_',num2str(ybgn(jj))) '/readEfield.m'])
            cd(strcat(outDielecdir_yzr,'_',num2str(ybgn(jj))))   
            [yzDielecEtensor]=readEfield(length(position_z),kx(ii),w(ii));
            DielecEtensor_yz(:,:,:,:,jj,ii)=yzDielecEtensor(:,:,:,:);
            cd(strcat(outVacdir_yzr,'_',num2str(ybgn(jj))))    
            [yzVacEtensor]=readEfield(length(position_z),kx(ii),w(ii));
            VacEtensor_yz(:,:,:,:,jj,ii)=yzVacEtensor(:,:,:,:);
        end
    end
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
                for l=1:size(DielecEtensor_yz,4)
                    for m=1:size(DielecEtensor_yz,5)
                        DielecEtensor_yzr(i,j,k,l,m)=2*sum(DielecEtensor_yz(i,j,k,l,m,:).*(mu0*c^2/a_meep^3));
                        VacEtensor_yzr(i,j,k,l,m)=2*sum(VacEtensor_yz(i,j,k,l,m,:).*(mu0*c^2/a_meep^3));
                    end
                end;
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
                 %Gamma_xy(i,j,l,kk)=sum(imag(g(:)).*real(DielecEtensor_xyr(:,i,j,l)).*dt.*a_meep/c);
                 %Gamma_vac_xy(i,j,l,kk)=sum(imag(g(:)).*real(VacEtensor_xyr(:,i,j,l)).*dt.*a_meep/c);
            end
            for p=1:length(position_z)
                for q=1:length(ybgn)
                    Gamma_yz(i,j,p,q,kk)=sum(imag(g(:)).*real(DielecEtensor_yzr(:,i,j,p,q)).*dt.*a_meep/c);
                    Gamma_vac_yz(i,j,p,q,kk)=sum(imag(g(:)).*real(VacEtensor_yzr(:,i,j,p,q)).*dt.*a_meep/c);
                end
            end
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
                E_shift_F_3_CP_xy(:,:,i,j,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                E_shift_F_3_Vac_xy(:,:,i,j,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);                               
                %E_shift_F_3_CP_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                %E_shift_F_3_Vac_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
                elseif ind_ff==2
                states2=[repmat(Fi(ind_ff),numel(mf_i1),1),mf_i1];
                [scalar_coeff_1,vector_coeff_1,tensor_coeff_1]=vectortensor_Casimir(states2,Fi(2),Green_xyr);
                [scalar_coeff_2,vector_coeff_2,tensor_coeff_2]=vectortensor_Casimir(states2,Fi(2),Vacuum_xyr);
                E_shift_F_4_CP_xy(:,:,i,j,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                E_shift_F_4_Vac_xy(:,:,i,j,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
                %E_shift_F_4_CP_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                %E_shift_F_4_Vac_xy(:,:,i,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
                end
            end
        end
        for l=1:length(position_z)
            for m=1:length(ybgn)
                Green_yzr=Gamma_yz(:,:,l,m,kk);
                Vacuum_yzr=Gamma_vac_yz(:,:,l,m,kk);
                if ind_ff==1
                    states1=[repmat(Fi(ind_ff),numel(mf_i1),1),mf_i1];
                    [scalar_coeff_1,vector_coeff_1,tensor_coeff_1]=vectortensor_Casimir(states1,Fi(1),Green_yzr);
                    [scalar_coeff_2,vector_coeff_2,tensor_coeff_2]=vectortensor_Casimir(states1,Fi(1),Vacuum_yzr);
                    E_shift_F_3_CP_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                    E_shift_F_3_Vac_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
                elseif ind_ff==2
                    states2=[repmat(Fi(ind_ff),numel(mf_i1),1),mf_i1];
                    [scalar_coeff_1,vector_coeff_1,tensor_coeff_1]=vectortensor_Casimir(states2,Fi(2),Green_yzr);
                    [scalar_coeff_2,vector_coeff_2,tensor_coeff_2]=vectortensor_Casimir(states2,Fi(2),Vacuum_yzr);
                    E_shift_F_4_CP_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_1-alpha1(kk,ind_ff).*vector_coeff_1-alpha2(kk,ind_ff).*tensor_coeff_1);
                    E_shift_F_4_Vac_yz(:,:,l,m,kk)=(-alpha0(kk,ind_ff).*scalar_coeff_2-alpha1(kk,ind_ff).*vector_coeff_2-alpha2(kk,ind_ff).*tensor_coeff_2);
                end
            end
        end
    end
end

%clear DielecEtensor_xyr DielecEtensor_zr VacEtensor_xyr VacEtensor_zr
%clear Gamma_xy Gamma_vac_xy Gamma_z Gamma_vac_z Green_xyr Green_zr Vacuum_xyr Vacuum_zr

%% Optical Stark shift from trap mode by MPB field profile
cd('/home/uqml/Desktop/MEEP/20170926_CP_final/Final_figure')
amplitude=7*10^6;
nband=22;
kpoint=2;
[HHXY_stark,HHYZ_stark]=Starkshift(xbgn,ybgn,position_y,position_z,amplitude,nband,kpoint,a_meep)
H_stark_xy_F3=HHXY_stark{1,1};
H_stark_xy_F4=HHXY_stark{1,2};
H_stark_yz_F3=HHYZ_stark{1,1};
H_stark_yz_F4=HHYZ_stark{1,2};
cd(workspace)
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
    
    [QTxy3,Total_xy_F3(:,:,k,l)]=eig(E_total_shift_F_3_CP_xy(:,:,k,l)+H_stark_xy_F3{k,l});
    [QTxy4,Total_xy_F4(:,:,k,l)]=eig(E_total_shift_F_4_CP_xy(:,:,k,l)+H_stark_xy_F4{k,l});
    Total_Eshift_xy(k,l,1)=trace(Total_xy_F3(:,:,k,l))-trace(Vac_xy_F3(:,:,k,l));
    Total_Eshift_xy(k,l,2)=trace(Total_xy_F4(:,:,k,l))-trace(Vac_xy_F4(:,:,k,l));
    Total_Potential_xy=Total_Eshift_xy(:,:,1)+Total_Eshift_xy(:,:,2);
    
       %for i=1:7
       %     for j=1:7
       %         E_total_shift_F_3_CP_xy(i,j,l)=sum(E_shift_F_3_CP_xy(i,j,l,:));
       %         E_total_shift_F_3_Vac_xy(i,j,l)=sum(E_shift_F_3_Vac_xy(i,j,l,:));
       %     end
       % end
       % for i=1:9
       %     for j=1:9
       %         E_total_shift_F_4_CP_xy(i,j,l)=sum(E_shift_F_4_CP_xy(i,j,l,:));
       %         E_total_shift_F_4_Vac_xy(i,j,l)=sum(E_shift_F_4_Vac_xy(i,j,l,:));
       %     end
       % end
    %[QCxy3,CP_xy_F3(:,:,l)]=eig(E_total_shift_F_3_CP_xy(:,:,l));
    %[QVxy3,Vac_xy_F3(:,:,l)]=eig(E_total_shift_F_3_Vac_xy(:,:,l));
    %[QCxy4,CP_xy_F4(:,:,l)]=eig(E_total_shift_F_4_CP_xy(:,:,l));
    %[QVxy4,Vac_xy_F4(:,:,l)]=eig(E_total_shift_F_4_Vac_xy(:,:,l));
    %CP_potential_xy(l,1)=trace(CP_xy_F3(:,:,l))-trace(Vac_xy_F3(:,:,l));
    %CP_potential_xy(l,2)=trace(CP_xy_F4(:,:,l))-trace(Vac_xy_F4(:,:,l));
    %Casimir_potential_xy=CP_potential_xy(:,1)+CP_potential_xy(:,2);
    end
end
save('CP_potential_xy.mat','Casimir_potential_xy')
%clear E_shift_F_3_CP_xy E_shift_F_3_Vac_xy E_shift_F_4_CP_xy E_shift_F_4_Vac_xy

for ii=1:length(position_z)
    for k=1:length(ybgn)
       for i=1:7
            for j=1:7
                E_total_shift_F_3_CP_yz(i,j,k,ii)=sum(E_shift_F_3_CP_yz(i,j,ii,k,:));
                E_total_shift_F_3_Vac_yz(i,j,k,ii)=sum(E_shift_F_3_Vac_yz(i,j,ii,k,:));
            end
        end
        for i=1:9
            for j=1:9
                E_total_shift_F_4_CP_yz(i,j,k,ii)=sum(E_shift_F_4_CP_yz(i,j,ii,k,:));
                E_total_shift_F_4_Vac_yz(i,j,k,ii)=sum(E_shift_F_4_Vac_yz(i,j,ii,k,:));
            end
        end    
    [QCyz3,CP_yz_F3(:,:,k,ii)]=eig(E_total_shift_F_3_CP_yz(:,:,k,ii));
    [QVyz3,Vac_yz_F3(:,:,k,ii)]=eig(E_total_shift_F_3_Vac_yz(:,:,k,ii));
    [QCyz4,CP_yz_F4(:,:,k,ii)]=eig(E_total_shift_F_4_CP_yz(:,:,k,ii));
    [QVyz4,Vac_yz_F4(:,:,k,ii)]=eig(E_total_shift_F_4_Vac_yz(:,:,k,ii));
    CP_potential_yz(k,ii,1)=trace(CP_yz_F3(:,:,k,ii))-trace(Vac_yz_F3(:,:,k,ii));
    CP_potential_yz(k,ii,2)=trace(CP_yz_F4(:,:,k,ii))-trace(Vac_yz_F4(:,:,k,ii));
    Casimir_potential_yz=CP_potential_yz(:,:,1)+CP_potential_yz(:,:,2);
    
    [QTyz3,Total_yz_F3(:,:,k,ii)]=eig(E_total_shift_F_3_CP_yz(:,:,k,ii)+H_stark_yz_F3{ii,k});
    [QTyz4,Total_yz_F4(:,:,k,ii)]=eig(E_total_shift_F_4_CP_yz(:,:,k,ii)+H_stark_yz_F4{ii,k});
    Total_Eshift_yz(k,ii,1)=trace(Total_yz_F3(:,:,k,ii))-trace(Vac_yz_F3(:,:,k,ii));
    Total_Eshift_yz(k,ii,2)=trace(Total_yz_F4(:,:,k,ii))-trace(Vac_yz_F4(:,:,k,ii));
    Total_Potential_yz=Total_Eshift_yz(:,:,1)+Total_Eshift_yz(:,:,2);
    
    end
end
save('CP_potential_yz.mat','Casimir_potential_yz')
%clear E_shift_F_3_CP_z E_shift_F_3_Vac_z E_shift_F_4_CP_z E_shift_F_4_Vac_z

%% Unit conversion // Matrix construction%%
CP_temp_unit_yz=Casimir_potential_yz./kB.*10^6/(6*pi);
CP_temp_unit_xy=Casimir_potential_xy./kB.*10^6/(6*pi);
Total_temp_unit_yz=Total_Potential_yz./kB.*10^6/(6*pi);
Total_temp_unit_xy=Total_Potential_xy./kB.*10^6/(6*pi);
%CP_freq_unit_xy=real(Casimir_potential_xy)./h.*10.^(-3) ; %units in KHz
CPZ=[fliplr(real(CP_temp_unit_yz(1,2:end))),real(CP_temp_unit_yz(1,:))];
TotalZ=[fliplr(real(Total_temp_unit_yz(1,2:end))),real(Total_temp_unit_yz(1,:))];
CPY=[fliplr(real(CP_temp_unit_xy(1,2:end))),real(CP_temp_unit_xy(1,:))];
TotalY=[fliplr(real(Total_temp_unit_xy(1,2:end))),real(Total_temp_unit_xy(1,:))];
CPX=[fliplr(real(CP_temp_unit_xy(2:end,1)')),real(CP_temp_unit_xy(:,1)')];
TotalX=[fliplr(real(Total_temp_unit_xy(2:end,1)')),real(Total_temp_unit_xy(:,1)')];
pxx=[-1*fliplr(xbgn(2:end))*a_meep*10^6,xbgn*a_meep*10^6];
pyy=[-1*fliplr(ybgn(2:end))*a_meep*10^6,ybgn*a_meep*10^6];
pzz=[-1*fliplr(position_z(2:end))*a_meep*10^6,position_z*a_meep*10^6];
CPYZ=[rot90(real(CP_temp_unit_yz(2:end,2:end)),2),flipud(real(CP_temp_unit_yz(2:end,:))); ...
    fliplr(real(CP_temp_unit_yz(:,2:end))),real(CP_temp_unit_yz)];
TotalYZ=[rot90(real(Total_temp_unit_yz(2:end,2:end)),2),flipud(real(Total_temp_unit_yz(2:end,:))); ...
    fliplr(real(Total_temp_unit_yz(:,2:end))),real(Total_temp_unit_yz)];
CPXY=[rot90(real(CP_temp_unit_xy(2:end,2:end)),2),flipud(real(CP_temp_unit_xy(2:end,:))); ...
    fliplr(real(CP_temp_unit_xy(:,2:end))),real(CP_temp_unit_xy)];
TotalXY=[rot90(real(Total_temp_unit_xy(2:end,2:end)),2),flipud(real(Total_temp_unit_xy(2:end,:))); ...
    fliplr(real(Total_temp_unit_xy(:,2:end))),real(Total_temp_unit_xy)];

%% Interporlation
Inpxx=linspace(-xbgn(end)*a_meep*10^6,xbgn(end)*a_meep*10^6,1000);
Inpyy=linspace(-ybgn(end)*a_meep*10^6,ybgn(end)*a_meep*10^6,1000);
Inpzz=linspace(-position_z(end)*a_meep*10^6,position_z(end)*a_meep*10^6,1000);

InCPX=interp1(pxx,CPX,Inpxx,'spline');
InCPY=interp1(pyy,CPY,Inpyy,'spline');
InCPZ=interp1(pzz,CPZ,Inpzz,'spline');
InTotalX=interp1(pxx,TotalX,Inpxx,'spline');
InTotalY=interp1(pyy,TotalY,Inpyy,'spline');
InTotalZ=interp1(pzz,TotalZ,Inpzz,'spline');

% [Xplotxy,Yplotxy]=meshgrid(pxx,pyy);
% [InXplotxy,InYplotxy]=meshgrid(Inpxx,Inpyy);
% [Yplotyz,Zplotyz]=meshgrid(pyy,pzz);
% [InYplotyz,InZplotyz]=meshgrid(Inpyy,Inpzz);
% InTotalXY=interp2(Xplotxy,Yplotxy,TotalXY',InXplotxy,InYplotxy,'linear');
% InTotalYZ=interp2(Yplotyz,Zplotyz,TotalYZ',InYplotyz,InYplotyz,'linear');

%% Colormap

    Mycolor=[0.15 0.15 0.85
        0.2 0.2 0.8
        0.2 0.4 0.8
        0.4 0.5 0.8
        0.8 0.5 0.4
        0.8 0.4 0.2
        0.8 0.2 0.2
        0.85 0.15 0.15]

    Blue=[0, 0, 0
        0, 0, 0.5
        0, 0, 1
        0, 0.5, 1
        0, 1, 1
        0.5,1,1
        1,0.5,0.5];
    
    Red=[0, 0, 0
        0.5, 0, 0
        1, 0, 0
        1, 0.5, 0
        1, 1, 0
        1,1,0.5
        1,1,1];

[X,Y]=meshgrid(1:size(Blue,1),1:size(Blue,2));
[X1,Y1]=meshgrid(1:0.01:size(Blue,1),1:size(Blue,2));
CBlue=interp2(X,Y,Blue',X1,Y1);
CRed=interp2(X,Y,Red',X1,Y1);
[X,Y]=meshgrid(1:size(Mycolor,1),1:size(Mycolor,2));
[X1,Y1]=meshgrid(1:0.01:size(Mycolor,1),1:size(Mycolor,2));
CMycolor=interp2(X,Y,Mycolor',X1,Y1);

%% Plot

figure(1)
h=surf(pzz,pyy,CPYZ)
hXLabel = xlabel(['z $(\mu m)$'],'position',[0 0.07 -120],'Rotation',8,'Fontsize',15,'Interpreter','Latex');
hYLabel = ylabel(['y $(\mu m)$'],'position',[0.48 -0.02 -120],'Rotation',-70,'Fontsize',15,'Interpreter','Latex');
hzLabel = zlabel(['$U_{\rm CP}$ $(\mu K)$'],'Fontsize',15,'Interpreter','Latex')
zlim([-80,0])
caxis([-80,0])
H=colorbar
set(H,'Position',[0.88 0.1100 0.0468 0.8150])
colormap(CMycolor')
set(h,'EdgeColor',[0 0 0])
shading interp
set(gca,'Fontsize',12)
set(gca,'ZTick',[-60:20:0])
set(gca,'YTick',[-0.1:0.05:0.1])
set(gca,'XTick',[-0.4:0.2:0.4])
set(gca,'Position',[0.15 0.15 0.65 0.75])
campos([90 70 40000])
grid on
save('CPYZ.mat','CPYZ')
save('pzz.mat','pzz')
save('pyy.mat','pyy')
saveas(h,['CPYZ.png'])

figure(2)
h=surf(pzz,pyy,TotalYZ)
hXLabel = xlabel(['z $(\mu m)$'],'position',[0.05 0.125 -120],'Rotation',14,'Fontsize',15,'Interpreter','Latex');
hYLabel = ylabel(['y $(\mu m)$'],'position',[0.55 0.03 -120],'Rotation',-68,'Fontsize',15,'Interpreter','Latex');
hzLabel = zlabel(['$U_{\rm total}$ $(\mu K)$'],'Fontsize',15,'Interpreter','Latex')
zlim([-100 20])
H=colorbar
set(H,'Position',[0.88 0.1100 0.0468 0.8150])
colormap(CMycolor')
caxis([-100 20])
set(h,'EdgeColor','none')
shading interp
set(gca,'Fontsize',12)
set(gca,'ZTick',[-100:100:100])
set(gca,'YTick',[-0.1:0.05:0.1])
set(gca,'XTick',[-0.4:0.2:0.4])
set(gca,'Position',[0.15 0.15 0.65 0.75])
campos([80 60 100000])
grid on
save('TOTALYZ.mat','TotalYZ')
saveas(h,['TOTALYZ.png'])

figure(3)
h1=plot(Inpzz,InTotalZ,'k-','linewidth',1.5)
hXLabel = xlabel(['z $(\mu m)$'],'Fontsize',16,'Interpreter','Latex');
hYLabel = ylabel(['$U_{\rm total}$ $(\mu K)$ '],'Fontsize',16,'Interpreter','Latex');
grid on
set(gca,'Fontsize',12)
set(gca,'YTick',[-40:5:0])
set(gca,'XTick',[-0.4:0.1:0.4])
save('InTotalZ.mat','InTotalZ')
save('Inpzz.mat','Inpzz')
saveas(h1,['TOTALZ.png'])
% 
figure(4)
h2=plot(Inpzz,InCPZ,'k-','linewidth',1.5)
hXLabel = xlabel(['z $(\mu m)$'],'Fontsize',16,'Interpreter','Latex');
hYLabel = ylabel(['$U_{\rm CP}$ $(\mu K)$ '],'Fontsize',16,'Interpreter','Latex');
grid on
set(gca,'Fontsize',12)
set(gca,'YTick',[-40:5:0])
set(gca,'XTick',[-0.4:0.1:0.4])
save('InCPZ.mat','InCPZ')
save('Inpzz.mat','Inpzz')
saveas(h2,['CPZ.png'])

figure(5)
h3=plot(Inpyy,InTotalY,'k-','linewidth',1.5)
hXLabel = xlabel(['y $(\mu m)$'],'Fontsize',16,'Interpreter','Latex');
hYLabel = ylabel(['$U_{\rm total}$ $(\mu K)$ '],'Fontsize',16,'Interpreter','Latex');
ylim([-100,40])
xlim([-0.08 0.08])
set(gca,'YTick',[-80:20:40])
set(gca,'XTick',[-0.08:0.04:0.08])
set(gca,'Fontsize',12)
grid on
save('InTotalY.mat','InTotalY')
save('Inpyy.mat','Inpyy')
saveas(h3,['TOTALY.png'])

figure(6)
h4=plot(Inpyy,InCPY,'k-','linewidth',1.5)
hXLabel = xlabel(['y $(\mu m)$'],'Fontsize',16,'Interpreter','Latex');
hYLabel = ylabel(['$U_{\rm CP}$ $(\mu K)$ '],'Fontsize',16,'Interpreter','Latex');
ylim([-100,0])
xlim([-0.08 0.08])
set(gca,'YTick',[-80:20:0])
set(gca,'XTick',[-0.08:0.04:0.08])
set(gca,'Fontsize',12)
grid on
save('InCPY.mat','InCPY')
save('Inpyy.mat','Inpyy')
saveas(h4,['CPY.png'])

figure(7)
h5=plot(Inpxx,InTotalX/1000,'k-','linewidth',1.5)
hXLabel = xlabel(['z $(\mu m)$'],'Fontsize',16,'Interpreter','Latex');
hYLabel = ylabel(['$U_{\rm total}$ $(mK)$ '],'Fontsize',16,'Interpreter','Latex');
set(gca,'Fontsize',12)
set(gca,'YTick',[0:4:16])
grid on
save('InTotalX.mat','InTotalX')
save('Inpxx.mat','Inpxx')
saveas(h5,['TOTALX.png'])

figure(8)
h6=plot(Inpxx,InCPX,'k-','linewidth',1.5)
hXLabel = xlabel(['z $(\mu m)$'],'Fontsize',16,'Interpreter','Latex');
hYLabel = ylabel(['$U_{\rm CP}$ $(\mu K)$ '],'Fontsize',16,'Interpreter','Latex');
ylim([[-44.1,-43.5]])
set(gca,'Fontsize',12)
set(gca,'YTick',[-44:0.1:-43.5])
grid on
save('Inpxx.mat','Inpxx')
save('InCPX.mat','InCPX')
saveas(h6,['CPX.png'])

figure(9)
h4=surf(pyy,pxx,CPXY/1000)
hXLabel = xlabel(['y $(\mu m)$'],'position',[0.04 0.37 -0.8],'Rotation',9,'Fontsize',15,'Interpreter','Latex');
hYLabel = ylabel(['x $(\mu m)$'],'position',[0.16 0.2 -0.8],'Rotation',-55,'Fontsize',15,'Interpreter','Latex');
hzLabel = zlabel(['$U_{\rm CP}$ (mK)'],'Fontsize',15,'Interpreter','Latex')
zlim([-1,0])
caxis([-1,0])
H=colorbar
set(H,'Position',[0.8900 0.1100 0.0468 0.8150])
colormap(CMycolor')
set(h4,'EdgeColor','none')
shading interp
set(gca,'Fontsize',12)
set(gca,'ZTick',[-0.8:0.4:0])
set(gca,'YTick',[-0.2:0.1:0.2])
set(gca,'XTick',[-0.1:0.05:0.1])
set(gca,'Position',[0.15 0.15 0.65 0.75])
campos([10 45 130])
grid on
save('CPXY.mat','CPXY')
save('pxx.mat','pxx')
save('pyy.mat','pyy')
saveas(h4,['CPXY.png'])

figure(10)
h5=surf(pyy,pxx,TotalXY/1000)
hXLabel = xlabel(['y $(\mu m)$'],'position',[0.015 0.25 -15],'Rotation',55,'Fontsize',15,'Interpreter','Latex');
hYLabel = ylabel(['x $(\mu m)$'],'position',[0.12 0.05 -15],'Rotation',-15,'Fontsize',15,'Interpreter','Latex');
hzLabel = zlabel(['$U_{\rm total}$ ($mK$)'],'Fontsize',15,'Interpreter','Latex')
H=colorbar
set(H,'Position',[0.8900 0.1100 0.0468 0.8150])
colormap(CMycolor')
caxis([-10 15])
set(h5,'EdgeColor','none')
shading interp
set(gca,'Fontsize',12)
set(gca,'ZTick',[-10:10:20])
set(gca,'Position',[0.15 0.15 0.65 0.75])
zlim([-10 20])
set(gca,'YTick',[-0.2:0.1:0.2])
set(gca,'XTick',[-0.1:0.05:0.1])
campos([15 15 3000])
grid on
save('TOTALXY.mat','TotalXY')
saveas(h5,['TOTALXY.png'])

%xx=position_y.*(a_meep*10.^6);
%figure(3)
%plot(xx,CP_freq_unit_xy,'bo')
%save('Sigma1.mat','CP_freq_unit_xy');
%hold on
%fun=@(c,x)(-1*c(2))./((x-c(3)).*x.^(3));
%c0=[0 0 0];
%c=lsqcurvefit(fun,c0,xx,CP_freq_unit_xy');
%x=linspace(xx(1),xx(end),1000);
%plot(x,fun(c,x),'r-','linewidth',1.5)
%hold off
%legend('Numercal data','curve fitting')
%xlabel('d (um)')
%ylabel('CP potential (kHz)')

%f=c(3) 
%C3=c(2)/c(3)    %1.16 kHz*um^3
%C4=c(2) %0.15 kHz*um^4

%figure(4)
%x1=linspace(0.001,100,100000);
%CP_planar=real(CP_freq_unit_xy').*xx.^3*1;
%func=@(C,X)(-1*C(2))./(X-C(3));
%C0=[0 0 0];
%C=lsqcurvefit(func,C0,xx,CP_planar);
%loglog(xx,CP_planar,'ro')
%hold on
%loglog(x1,-1*C(2)./(x1-C(3)),'k-')
%%hold off
%legend('Numercal data','curve fitting')
%xlabel('d (um)')
%ylabel('Ucp*d^3 (kHz*um^3)')

% CC3=C(2)/C(3)
% CC4=C(2)
