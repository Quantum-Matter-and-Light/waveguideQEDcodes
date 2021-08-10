%%

function [ffield,Gx,Gy,Gz]=Efieldfrequencyplot(GFresult,zeropos,mu)
cc=3e8;eps0=8.854e-12;


Efield=GFresult.E;
ffield=GFresult.f;
w=2*pi*ffield;
zeroindex=zeropos;
%round(zeropos(1)*zeropos(2));
Gx=zeros(size(ffield));
Gy=zeros(size(ffield));
Gz=zeros(size(ffield));

for i1=1:size(ffield)
        Gx(i1)=squeeze((Efield(zeroindex,1,i1))).*cc.^2.*eps0./w(i1).^2./mu;
        Gy(i1)=squeeze((Efield(zeroindex,2,i1))).*cc.^2.*eps0./w(i1).^2./mu;
        Gz(i1)=squeeze((Efield(zeroindex,3,i1))).*cc.^2.*eps0./w(i1).^2./mu;
end

