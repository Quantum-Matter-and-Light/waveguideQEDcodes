%%
function planar(h,cont,nSiN,tSiN)
lcommand(h,'addrect',cont);
lset(h,cont,'name',['myrect']);
lset(h,cont,'x min',-0.4e-6);
lset(h,cont,'x max',0);
lset(h,cont,'y min',-1.0e-6);
lset(h,cont,'y max',1.0e-6);
lset(h,cont,'z min',-1.0e-6);
lset(h,cont,'z max',1.0e-6)'
lset(h,cont,'index',nSiN);
%lset(h,cont,'x min bc','pml');
%lset(h,cont,'x max bc','pml');
end