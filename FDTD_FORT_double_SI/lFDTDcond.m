%%
function lFDTDcond(h,cont,meshacc,timefactor,SimX,SimY,SimZ,xSimMin,xSimMax,ySimMin,ySimMax,zSimMin,zSimMax,dx,dy,dz,AutoShutOff,Courant)
cc=299792458; %speed of light
Tf=timefactor*10e-6/cc;

lcommand(h,'addfdtd',cont)
%lset(h,cont,'dimension',1)  % 1 = 2D, 2= 3D
lset(h,cont,'mesh accuracy',meshacc);
lset(h,cont,'simulation time',Tf);
lset(h,cont,'x',SimX);
lset(h,cont,'y',SimY);
lset(h,cont,'z',SimZ);
lset(h,cont,'x min',xSimMin);lset(h,cont,'x max',xSimMax);
lset(h,cont,'y min',ySimMin);lset(h,cont,'y max',ySimMax);
lset(h,cont,'z min',zSimMin);lset(h,cont,'z max',zSimMax);
lset(h,cont,'auto shutoff min',AutoShutOff);
lset(h,cont,'mesh type','uniform');
lset(h,cont,'min mesh step',2.5e-12);
lset(h,cont,'define x mesh by','maximum mesh step');
lset(h,cont,'define y mesh by','maximum mesh step');
lset(h,cont,'define z mesh by','maximum mesh step');
lset(h,cont,'dx',dx);
lset(h,cont,'dy',dy);
lset(h,cont,'dz',dz);
lset(h,cont,'dt stability factor',Courant);
lset(h,cont,'mesh refinement','conformal variant 0');

delete('FDTD_Condition.txt');
fileID=fopen('FDTD_Condition.txt','w');
Parameters=['FDTD_Condition \r\n' ...
    'time factor: ' num2str(timefactor) '\r\n' ...
    'Simulation time: ' num2str(Tf) 's \r\n' ...
    'Mesh accuracy: ' num2str(meshacc) 'm \r\n' ...
    'FDTD_SimX: ' num2str(SimX) 'm \r\n'...
    'FDTD_SimY: ' num2str(SimY) 'm \r\n'...
    'FDTD_SimZ: ' num2str(SimZ) 'm \r\n'...
    'FDTD_SimXspan: ' num2str(2*xSimMax) 'm \r\n'...
    'FDTD_SimYspan: ' num2str(2*ySimMax) 'm \r\n'...
    'FDTD_SimZspan: ' num2str(2*zSimMax) 'm \r\n'...
    'AutoShutOff_minimum: ' num2str(AutoShutOff) '\r\n' ...
    'Meshtype: uniform' '\r\n' ...
    'Resolution (dx=dy=dz): ' num2str(dx) '\r\n' ...
    'Time Stability Factor (Courant): ' num2str(Courant) '\r\n' ...
    'Symmetric BC?:' 'Yes' '\r\n'];
fprintf(fileID,Parameters);
fclose(fileID);

end
