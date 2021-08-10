function lSTLimport(h,cont,PATH,nSiN,x,y,z,toggle)
if toggle==1;
    lcommand(h,['stlimport("' PATH '")'],cont);
    lset(h,cont,'name','SPCW_structure');
    lset(h,cont,'index',nSiN);    
    lset(h,cont,'x',x);
    lset(h,cont,'y',y);
    lset(h,cont,'z',z);
elseif toggle==2;
    lcommand(h,'select("SPCW_structure")',true);
    lcommand(h,'delete',true);
end
end