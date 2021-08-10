%%
function lFDTDmonitors(h,cont,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,xpos,ypos,zpos,freqhighres,timefactor)
cc=299792458; %speed of light
Tf=timefactor*10e-6/cc;

lcommand(h,'addpower',cont)
lset(h,cont,'name','GreenFunctionDATA');

lset(h,cont,'monitor type','linear x');
%     lset(h,cont,'monitor type','2D Z-normal');
lset(h,cont,'x min',xSimMin);lset(h,cont,'x max',xSimMax);
%     lset(h,cont,'y min',ySimMin);lset(h,cont,'y max',ySimMax);
%     lset(h,cont,'x',xpos);
lset(h,cont,'y',ypos);
lset(h,cont,'z',zpos);

lset(h,cont,'override global monitor settings',1);
lset(h,cont,'frequency points',freqhighres);
lset(h,cont,'spatial interpolation','none');
lset(h,cont,'apodization','end');
% lset(h,cont,'apodization center',Tf/2);
lset(h,cont,'apodization time width',Tf/2);

end