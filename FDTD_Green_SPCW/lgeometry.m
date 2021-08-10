function lgeometry(h,nSiN,dipoleaxis,dx,dy,dz,toggle)
%% Script to make the W1PCW structure
tSiN=200e-9;
wid=228e-9;
thk=200e-9;
a0=370e-9;
ny=8;
r1=100e-9;
r2=106e-9;
r3=110e-9;
s1=-10e-9;
s2=0;
offset=105e-9;
Alpha=1.256;

delete('Parameters.txt');
fileID=fopen('Parameters.txt','w');
Parameters=['Structure Parameters \r\n' ...
    'Lattice Constant: ' num2str(a0) 'm \r\n' ...
    'Slab Thickness: ' num2str(thk) 'm \r\n' ...
    'Slot Width: ' num2str(wid) 'm \r\n' ...
    'r1 (First hole radius): ' num2str(r1) 'm \r\n'...
    'r2 (Second hole radius): ' num2str(r2) 'm \r\n'...
    'r3 (All the other hole radius): ' num2str(r3) 'm \r\n'...
    'S1 (Center-shift of the First hole): ' num2str(s1) 'm \r\n'...
    'S2 (Center-shift of the Second hole): ' num2str(s2) 'm \r\n'...
    'Offset (Add dielectric at the edge of the slot): ' num2str(offset) 'm \r\n'...
    'Alpha (Squircle): ' num2str(Alpha) '\r\n'];
fprintf(fileID,Parameters);
fclose(fileID);

if dipoleaxis=='x'
    x0=a0/2+dx/2;
    y0=0;
    z0=0;
elseif dipoleaxis=='y'
    x0=a0/2;
    y0=dy/2;
    z0=0;
elseif dipoleaxis=='z'
    x0=a0/2;
    y0=0;
    z0=dz/2;
end

if toggle==1;
    W1PCWunitcell(h,true,nSiN,tSiN,x0,y0,z0,a0,wid,thk,ny,r1,r2,r3,s1,s2,offset,Alpha);
elseif toggle==2;
    lcommand(h,'select("W1PhCW_Structure")',true);
    lcommand(h,'delete',true);
end

%% homogeneous
%planar(h,true,nSiN,tSiN)
%cylinder(h,true,nSiN,tSiN)
end

