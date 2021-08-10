function [Fi,mfi,omega_ff,Fi_easy,p_index,FF,alpha1_coeff,alpha2_coeff]=CP_groundstate_shift()
ep0=8.85e-12;
hbar=1.05457148e-34;
c=299792458;
I=7/2;
Li=0; %S_1/2
[Li,Ji,Fi,mfi]=spdf(Li);
Lf=1; %P_3/2
[Lf,Jf,Ff,mff]=spdf(Lf);
Fi=Fi{1};

pshalf=importdata('pshalf.mat'); %wavelength and life time for D1 and D2 line
p_index=pshalf(:,1); %principle quantum number
omega_ps=(2*pi*c)./(pshalf(:,[2 4])*1e-9); %transition frequencies for D1 (first column) and D2 (second column) line
tau_ps=pshalf(:,[3,5])*1e-6; %lifetime for D1 and D2 line
gamma_ps=1./tau_ps; %spontaneous emission rates

jj_factor=repmat((2*Ji+1)./(2*Jf+1),size(pshalf,1),1);
jj_ps_sq=((3*pi*ep0*hbar*c^3)./(omega_ps.^3).*(1./(jj_factor)).*gamma_ps); 

Ff_eas=[Ff{1},Ff{2}];
Fi_easy=repmat(Fi',1,numel(Ff_eas));
Ff_easy=repmat(Ff_eas,numel(Fi),1);
Ji_easy=repmat([Ji;Ji],1,numel(Ff_eas));
Jf_easy=repmat([Jf(1)*ones(size(Ff{1})),Jf(2)*ones(size(Ff{2}))],numel(Jf),1);
FF=zeros(size(Ff_easy,1),size(Ff_easy,2));
omega_ff=zeros(size(FF));
alpha1_coeff=zeros(size(FF));
alpha2_coeff=zeros(size(FF));

for ind_j=1:numel(p_index)
    for ind_f=1:numel(Ff_eas)
        %for Fi=3
        if Jf_easy(1,ind_f)==Jf(1)
            JdJ=sqrt(jj_ps_sq(ind_j,1));
        else
            JdJ=sqrt(jj_ps_sq(ind_j,2));
        end
        FF(1,ind_f,ind_j)=JdJ...
            *(-1)^(Ff_easy(1,ind_f)+Ji_easy(1,ind_f)+1+I)...
            *sqrt((2*Ff_easy(1,ind_f)+1)*(2*Ji_easy(1,ind_f)+1))...
            *Wigner6jcoeff(Ji_easy(1,ind_f),Jf_easy(1,ind_f),1,Ff_easy(1,ind_f),Fi_easy(1,ind_f),I);
        
        alpha1_coeff(1,ind_f,ind_j)=(-1)^(Fi_easy(1,ind_f)+Ff_easy(1,ind_f)+1)...
            *sqrt(6*Fi_easy(1,ind_f)*(2*Fi_easy(1,ind_f)+1)/(Fi_easy(1,ind_f)+1))...
            *Wigner6jcoeff(1,1,1,Fi_easy(1,ind_f),Fi_easy(1,ind_f),Ff_easy(1,ind_f));
        
        alpha2_coeff(1,ind_f,ind_j)=(-1)^(Ff_easy(1,ind_f)+Fi_easy(1,ind_f))...
            *sqrt(40*Fi_easy(1,ind_f)*(2*Fi_easy(1,ind_f)+1)...
            *(2*Fi_easy(1,ind_f)-1)/(3*(Fi_easy(1,ind_f)+1)*(2*Fi_easy(1,ind_f)+3)))...
            *Wigner6jcoeff(1,1,2,Fi_easy(1,ind_f),Fi_easy(1,ind_f),Ff_easy(1,ind_f));
        
        % for Fi=4
        FF(2,ind_f,ind_j)=JdJ...
            *(-1)^(Ff_easy(2,ind_f)+Ji_easy(2,ind_f)+1+I)...
            *sqrt((2*Ff_easy(2,ind_f)+1)*(2*Ji_easy(2,ind_f)+1))...
            *Wigner6jcoeff(Ji_easy(2,ind_f),Jf_easy(2,ind_f),1,Ff_easy(2,ind_f),Fi_easy(2,ind_f),I);
        
        alpha1_coeff(2,ind_f,ind_j)=(-1)^(Fi_easy(2,ind_f)+Ff_easy(2,ind_f)+1)...
            *sqrt(6*Fi_easy(2,ind_f)*(2*Fi_easy(2,ind_f)+1)/(Fi_easy(2,ind_f)+1))...
            *Wigner6jcoeff(1,1,1,Fi_easy(2,ind_f),Fi_easy(2,ind_f),Ff_easy(2,ind_f));
        
        alpha2_coeff(2,ind_f,ind_j)=(-1)^(Ff_easy(2,ind_f)+Fi_easy(2,ind_f))...
            *sqrt(40*Fi_easy(2,ind_f)*(2*Fi_easy(2,ind_f)+1)...
            *(2*Fi_easy(2,ind_f)-1)/(3*(Fi_easy(2,ind_f)+1)*(2*Fi_easy(2,ind_f)+3)))...
            *Wigner6jcoeff(1,1,2,Fi_easy(2,ind_f),Fi_easy(2,ind_f),Ff_easy(2,ind_f));
        
        % fill up the omega_ff matrix
        %for Fi=3
        if Jf_easy(1,ind_f)==Jf(1)
            omega_ff(1,ind_f,ind_j)=omega_ps(ind_j,1);
        elseif Jf_easy(1,ind_f)==Jf(2)
            omega_ff(1,ind_f,ind_j)=omega_ps(ind_j,2);
        end
        % for Fi=4
        if Jf_easy(2,ind_f)==Jf(1)
            omega_ff(2,ind_f,ind_j)=omega_ps(ind_j,1);
        elseif Jf_easy(2,ind_f)==Jf(2)
            omega_ff(2,ind_f,ind_j)=omega_ps(ind_j,2);
        end
    end
end
