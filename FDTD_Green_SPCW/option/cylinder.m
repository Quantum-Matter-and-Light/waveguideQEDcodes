function cylinder(h,cont,nSiN,tSiN)
lcommand(h,'addcircle',cont);
lset(h,cont,'name',['mycyl']);
lset(h,cont,'x',-0.2e-6);
lset(h,cont,'y',0);
lset(h,cont,'z',0);
lset(h,cont,'radius',0.2e-6);
lset(h,cont,'z span',3.2e-6);
lset(h,cont,'index',nSiN);
%lset(h,cont,'x min bc','pml');
%lset(h,cont,'x max bc','pml');
end