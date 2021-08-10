function [FdF]=DipoleZeeman(Jg,Je,Fg,Fe,me,mg,I,q)
%Jg=1/2;Je=3/2;I=7/2;Fe=5;Fg=4;me=1;mg=2;q=-1;
F=ClebschGordan(Fg,1,Fe,mg,q,me)/sqrt(2*Fe+1);
FdF=(-1)^(Fe+Jg+1+I)*sqrt((2*Fg+1)*(2*Je+1)) ...
        *Wigner6jcoeff(Fg,1,Fe,Je,I,Jg)*F/sqrt(2*Je+1);


