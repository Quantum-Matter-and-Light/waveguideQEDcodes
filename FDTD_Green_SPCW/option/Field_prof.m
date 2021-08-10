function [Xplot,Yplot,Gx,Gy,Gz,zeropos]=Field_prof(GFresult,mu,f,targetx,targetz)
cc=3e8;eps0=8.854e-12;
Efield=GFresult.E;
ffield=GFresult.f;
tmp=abs(ffield-f);
[idx idx]=min(tmp);
closest=ffield(idx);
w=2*pi*closest;
[idx1 idx1]=min(abs(GFresult.x-targetx));
[idx2 idx2]=min(abs(GFresult.z-targetz));
zeropos=idx1*idx2;

[Xplot,Yplot]=meshgrid(GFresult.x,GFresult.z);
[Xsize,Ysize]=size(Xplot);
Gx=zeros(size(Xplot));
irun=0;

%% Photonic crystal structure (2D measurements: monitor-type "2D y-normal") %%
for i1=1:Xsize
    for i2=1:Ysize
        irun=irun+1;
        Gx(i1,i2)=squeeze((Efield(irun,1,idx)))*cc^2*eps0/w^2/mu;
    end
end

Gy=zeros(size(Xplot));
irun=0;
for i1=1:Xsize
    for i2=1:Ysize
        irun=irun+1;
        Gy(i1,i2)=squeeze((Efield(irun,2,idx)))*cc^2*eps0/w^2/mu;
    end
end

Gz=zeros(size(Xplot));
irun=0;
for i1=1:Xsize
    for i2=1:Ysize
        irun=irun+1;
        Gz(i1,i2)=squeeze((Efield(irun,3,idx)))*cc^2*eps0/w^2/mu;
    end
end