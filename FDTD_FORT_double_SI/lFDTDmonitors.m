%%
function lFDTDmonitors(h,cont,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,xpos,ypos,zpos,freqhighres)
lcommand(h,'addpower',cont)
lset(h,cont,'name','SILDATA_X');
lset(h,cont,'monitor type','2D X-normal');
lset(h,cont,'y min',ySimMin);lset(h,cont,'y max',ySimMax);
lset(h,cont,'z min',zSimMin);lset(h,cont,'z max',zSimMax);
lset(h,cont,'x',xpos);
lset(h,cont,'y',ypos);
lset(h,cont,'z',zpos);
lset(h,cont,'override global monitor settings',1);
lset(h,cont,'frequency points',freqhighres);
lset(h,cont,'spatial interpolation','none');

lcommand(h,'addpower',cont)
lset(h,cont,'name','SILDATA_Z');
lset(h,cont,'monitor type','2D Z-normal');
lset(h,cont,'x min',xSimMin/8);lset(h,cont,'x max',xSimMax/8);
lset(h,cont,'y min',ySimMin);lset(h,cont,'y max',ySimMax);
lset(h,cont,'x',xpos);
lset(h,cont,'y',ypos);
lset(h,cont,'z',zpos);
lset(h,cont,'override global monitor settings',1);
lset(h,cont,'frequency points',freqhighres);
lset(h,cont,'spatial interpolation','none');

lcommand(h,'addpower',cont)
lset(h,cont,'name','SILDATA_Y');
lset(h,cont,'monitor type','2D Y-normal');
lset(h,cont,'x min',xSimMin);lset(h,cont,'x max',xSimMax);
lset(h,cont,'z min',zSimMin);lset(h,cont,'z max',zSimMax);
lset(h,cont,'x',xpos);
lset(h,cont,'y',ypos);
lset(h,cont,'z',zpos);
lset(h,cont,'override global monitor settings',1);
lset(h,cont,'frequency points',freqhighres);
lset(h,cont,'spatial interpolation','none');

end