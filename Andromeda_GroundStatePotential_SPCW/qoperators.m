function [coeff_f,state_f]=qoperators(F,state_i)
%F = +1 (F+), -1 (F-), 0 (Fz), 2 (F^2)
%state_i is |F mf>
%state_i=[3 -3];
%F=-1
hbar=1;
f=state_i(1);
mf=state_i(2);
if F==1
    coeff_f=hbar*sqrt(f*(f+1)-mf*(mf+1));
    if mf<f
    state_f=[f,mf+1];
    else
        state_f=state_i;
    end
elseif F==-1
    coeff_f=hbar*sqrt(f*(f+1)-mf*(mf-1));
    if -mf<f
    state_f=[f,mf-1];
        else
        state_f=state_i;
    end
elseif F==0
    coeff_f=hbar*mf;
    state_f=state_i;
elseif F==2
    coeff_f=hbar^2*f*(f+1);
    state_f=state_i;
else
    error('wrong format');
end
    
    
