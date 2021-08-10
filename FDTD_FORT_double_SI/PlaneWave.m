function []=PlaneWave(h,cont,name,freq,bw,x,y,z,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,axis,Atheta,Aphi,PolAngle,direction)
lcommand(h,'addplane',cont);
lset(h,cont,'name',name);
lset(h,cont,'injection axis',axis);
lset(h,cont,'direction',direction);
lset(h,cont,'angle theta',Atheta);
lset(h,cont,'angle phi',Aphi);
lset(h,cont,'polarization angle',PolAngle);
lset(h,cont,'x',x);lset(h,cont,'y',y);lset(h,cont,'z',z);
lset(h,cont,'x max',xSimMax);lset(h,cont,'x min',xSimMin);
lset(h,cont,'y max',ySimMax);lset(h,cont,'y min',ySimMin);
% lset(h,cont,'z max',zSimMax);lset(h,cont,'z min',zSimMin);
lset(h,cont,'center frequency',freq);
lset(h,cont,'frequency span',bw);
end