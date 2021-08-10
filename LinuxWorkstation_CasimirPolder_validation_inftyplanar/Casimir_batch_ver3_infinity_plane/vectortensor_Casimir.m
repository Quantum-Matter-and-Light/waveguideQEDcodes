%% Spherical tensor decomposition

function [scalar_coeff,vector_coeff,tensor_coeff]=vectortensor_Casimir(states,F,Green)
%tensor test
%Based on Solution from Mathematica code by Kyung and Aki
%Note that for vector_coeff, the 1i in front of the vector shift has
%already been taken into account here

%F=3;
n=2*F+1;

tensor_coeff=zeros(n);
vector_coeff=zeros(n);
scalar_coeff=zeros(n);

%states=[4  -4; 4  -3 ; 4  -2 ; 4  -1 ; 4  0 ; 4  1 ; 4  2 ; 4  3 ; 4  4];
%states=[3  -3; 3  -2 ; 3  -1 ; 3   0 ; 3  1 ; 3  2 ; 3  3];

%E=[1/2+1i,1+(1/2)*1i,1/3]

%Ex=(1/sqrt(2))*((E(1)/2)-(E(2)/2))
%Ey=(1i/sqrt(2))*((E(1)/2)+(E(2)/2))
%Ez=E(3)/2
%Green=[Ex*conj(Ex) Ex*conj(Ey) Ex*conj(Ez);Ey*conj(Ex) Ey*conj(Ey) Ey*conj(Ez);Ez*conj(Ex) Ez*conj(Ey) Ez*conj(Ez)]
%Green=[1/2 j/sqrt(2) -1/sqrt(2);-j/sqrt(2) 1/2 j/sqrt(2);-1/sqrt(2) -j/sqrt(2) 1]
%E=[1 0 0];
%Green=E'*E;
%G=importdata('G.mat')
%Green=G(:,:,1)
%Spherical tensor representation of Green tensor

G00=(1/4)*trace(Green);
G1p1=(1/2)*((Green(1,3)-Green(3,1))+1i*(Green(2,3)-Green(3,2)));
G1m1=(-1/2)*((Green(1,3)-Green(3,1))-1i*(Green(2,3)-Green(3,2)));
G1z0=-1i*(Green(1,2)-Green(2,1));
G2p2=(1/2)*((Green(1,1)-Green(2,2))+1i*(Green(1,2)+Green(2,1)));
G2m2=(1/2)*((Green(1,1)-Green(2,2))-1i*(Green(1,2)+Green(2,1)));
G2p1=(1/2)*((Green(3,1)+Green(1,3))+1i*(Green(3,2)+Green(2,3)));
G2m1=(1/2)*((Green(3,1)+Green(1,3))-1i*(Green(3,2)+Green(2,3)));
G2z1=(Green(1,1)+Green(2,2));
G2z2=Green(3,3);

scalar_coeff=G00*eye(n); %check!

for kk=1:n
    state_left=states(kk,:);
    for ll=1:n
        state_i=states(ll,:);
        %% Vector coefficient
        [coeff_fp,state_fp]=qoperators(1,state_i); %Fp terms 
        if state_left==state_fp
           vector_co_p=(G1p1)*coeff_fp/(4*F);
        else
           vector_co_p=0;
        end
        [coeff_fm,state_fm]=qoperators(-1,state_i); %Fm terms
        if state_left==state_fm
           vector_co_m=(G1m1)*coeff_fm/(4*F);
        else
           vector_co_m=0;
        end
        [coeff_fz,state_fz]=qoperators(0,state_i); %Fz terms
        if state_left==state_fz
           vector_co_z=(G1z0*coeff_fz)/(4*F);
        else
           vector_co_z=0;
        end
        vector_coeff(kk,ll)=vector_co_p+vector_co_m+vector_co_z;
	
	%% Tensor coefficient
        [coeff_m,state_m]=qoperators(-1,state_i); %Fm term
        if state_left==state_m
           tensor_co_m=((3*G2m1)/(8*F*(-1 + 2*F)))*coeff_m;
        else
           tensor_co_m=0;
        end
        
        [coeff_m2,state_m2]=qoperators(-1,state_m); %Fm^2
        if state_left==state_m2
           tensor_co_m2=((3*G2m2)/(8*F*(-1 + 2*F)))*coeff_m2*coeff_m;
        else
	       tensor_co_m2=0;
        end
        
        [coeff_p,state_p]=qoperators(1,state_i); %Fp term
        if state_left==state_p
           tensor_co_p=((-3*G2p1)/(8*F*(-1 + 2*F)))*coeff_p;
        else
           tensor_co_p=0;
        end
        
        [coeff_p2,state_p2]=qoperators(1,state_p); %Fp^2
        if state_left==state_p2
	       tensor_co_p2=((3*G2p2)/(8*F*(-1 + 2*F)))*coeff_p2*coeff_p;
        else
           tensor_co_p2=0;
        end
        
        [coeff_z,state_z]=qoperators(0,state_i); %Fz term
        if state_left==state_z
           tensor_co_z=(3*G2z1)/(8*F*(-1 + 2*F))*coeff_z;
        else
	       tensor_co_z=0;
        end
        
        [coeff_z2,state_z2]=qoperators(0,state_z); %Fz^2 term
        if state_left==state_z2
           tensor_co_z2=(3*G2z2)/(4*F*(-1 + 2*F))*coeff_z*coeff_z2;
        else
           tensor_co_z2=0;
        end

        if state_left==state_i
            tensor_diag=-abs(Green(1,1))/(4*(-1 + 2*F))-(F*(abs(Green(1,1))))/(4*(-1 + 2*F))...
             -(abs(Green(2,2)))/(4*(-1 + 2*F))-(F*(abs(Green(2,2))))/(4*(-1 + 2*F))...
             -(abs(Green(3,3)))/(4*(-1 + 2*F))-(F*(abs(Green(3,3))))/(4*(-1 + 2*F));
        else
           tensor_diag=0;
        end

        [coeff_pm_1,state_pm_1]=qoperators(1,state_i); %F+
        [coeff_pm_2,state_pm_2]=qoperators(-1,state_pm_1); %F-
        if state_left==state_pm_2
           tensor_co_pm=(3*G2z1)/(8*F*(-1 + 2*F))*coeff_pm_1*coeff_pm_2;
        else
           tensor_co_pm=0;
        end
        
        [coeff_mz_1,state_mz_1]=qoperators(0,state_i); %Fz
        [coeff_mz_2,state_mz_2]=qoperators(-1,state_mz_1); %F-
        if state_left==state_mz_2
           tensor_co_mz=((-3*G2m1)/(4*F*(-1+2*F)))*coeff_mz_1*coeff_mz_2;
        else
           tensor_co_mz=0;
        end

        [coeff_pz_1,state_pz_1]=qoperators(0,state_i); %Fz
        [coeff_pz_2,state_pz_2]=qoperators(1,state_pz_1); %F+
        if state_left==state_pz_2
           tensor_co_pz=((-3*G2p1)/(4*F*(-1+2*F)))*coeff_pz_1*coeff_pz_2;
        else
           tensor_co_pz=0;
        end
    
           tensor_coeff(kk,ll)=tensor_co_m+tensor_co_m2+tensor_co_p+tensor_co_p2+...
            tensor_co_z+tensor_co_z2+tensor_diag...
            +tensor_co_pm+tensor_co_mz+tensor_co_pz;
        end
    end
%tensor_coeff
