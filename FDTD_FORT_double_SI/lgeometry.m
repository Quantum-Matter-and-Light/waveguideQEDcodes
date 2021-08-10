function lgeometry(h,nSiN,dx,dy,dz,toggle)
%% Script to make the W1PCW structure
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
Ncell=10;

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
    'Alpha (Squircle): ' num2str(Alpha) '\r\n' ...
    'Number of Unit cell:' num2str(Ncell) '\r\n'];
fprintf(fileID,Parameters);
fclose(fileID);

if toggle==1;
    W1PCWunitcell(h,true,nSiN,0,0,0,a0,wid,thk,ny,r1,r2,r3,s1,s2,offset,Alpha,Ncell);
elseif toggle==2;
    lcommand(h,'selectpartial("W1PhCW_cell")',true);
    lcommand(h,'delete',true);
end

%% homogeneous
%planar(h,true,nSiN,tSiN)
%cylinder(h,true,nSiN,tSiN)
end

