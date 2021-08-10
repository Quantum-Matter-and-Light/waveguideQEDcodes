%%
function APCWunitcell(h,cont,nSiN,tSiN,w1,w2,a0,Ain,x0,gap,ncell)

lcommand(h,'addcustom',cont);
lset(h,cont,'name',['APCW section L ' num2str(round(x0/a0/ncell))]);
lset(h,cont,'create 3D object by','extrusion');
lset(h,cont,'equation units','m'); 
lset(h,cont,'index',nSiN);
lset(h,cont,'make nonsymmetric',1);
lset(h,cont,'equation 1',[num2str(w1) '+(' num2str(w2) '-' num2str(w1) ')*(x-' num2str(x0) ')/' num2str(a0) '+(' num2str(Ain) '/2)*' '(1+cos(2*pi*(x-' num2str(x0) ')/(' num2str(a0) ')))']);
lset(h,cont,'equation 2',0);
lset(h,cont,'x min',-a0/2+x0);
lset(h,cont,'x max',a0/2+a0*ncell+x0);
lset(h,cont,'y',gap/2);
lset(h,cont,'y span',1000e-9); % computational size
lset(h,cont,'z span',tSiN);
lset(h,cont,'equation 1',[num2str(w1) '+(' num2str(w2) '-' num2str(w1) ')*(x-' num2str(x0) ')/' num2str(a0) '+(' num2str(Ain) '/2)*' '(1+cos(2*pi*(x-' num2str(x0) ')/(' num2str(a0) ')))']);
lset(h,cont,'equation 2',0);
lset(h,cont,'detail',1)

lcommand(h,'addcustom',cont);
lset(h,cont,'name',['APCW section R ' num2str(round(x0/a0/ncell))]);
lset(h,cont,'create 3D object by','extrusion');
lset(h,cont,'equation units','m'); 
lset(h,cont,'index',nSiN);
lset(h,cont,'make nonsymmetric',1);
lset(h,cont,'equation 1',[num2str(w1) '+(' num2str(w2) '-' num2str(w1) ')*(x-' num2str(x0) ')/' num2str(a0) '+(' num2str(Ain) '/2)*' '(1+cos(2*pi*(x-' num2str(x0) ')/(' num2str(a0) ')))']);
lset(h,cont,'equation 2',0);
lset(h,cont,'x min',-a0/2+x0);
lset(h,cont,'x max',a0/2+a0*ncell+x0);
lset(h,cont,'y',-gap/2);
lset(h,cont,'y span',1000e-9); % computational size
lset(h,cont,'z span',tSiN);
lset(h,cont,'equation 1',[num2str(w1) '+(' num2str(w2) '-' num2str(w1) ')*(x-' num2str(x0) ')/' num2str(a0) '+(' num2str(Ain) '/2)*' '(1+cos(2*pi*(x-' num2str(x0) ')/(' num2str(a0) ')))']);
lset(h,cont,'equation 2',0);
lset(h,cont,'first axis','x');
lset(h,cont,'rotation 1',180);
lset(h,cont,'detail',1)
end