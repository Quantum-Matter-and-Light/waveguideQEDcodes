%%

function lindex(h,cont,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,xpos,ypos,zpos)

lcommand(h,'addindex',cont);
lset(h,cont,'name','index_monitor1');
lset(h,cont,'monitor type','2D X-normal');
lset(h,cont,'x',xpos);
% lset(h,cont,'x min',xSimMin);lset(h,cont,'x max',xSimMax);
lset(h,cont,'y',ypos);
lset(h,cont,'y min',ySimMin);lset(h,cont,'y max',ySimMax);
lset(h,cont,'z',zpos);
lset(h,cont,'z min',zSimMin);lset(h,cont,'z max',zSimMax);

lcommand(h,'addindex',cont);
lset(h,cont,'name','index_monitor2');
lset(h,cont,'monitor type','2D Z-normal');
lset(h,cont,'x',xpos);
lset(h,cont,'x min',xSimMin);lset(h,cont,'x max',xSimMax);
lset(h,cont,'y',ypos);
lset(h,cont,'y min',ySimMin);lset(h,cont,'y max',ySimMax);
lset(h,cont,'z',zpos);
% lset(h,cont,'z min',zSimMin);lset(h,cont,'z max',zSimMax);

end