function [Etensor]=readEfieldNk(Ny,kx)
w=1;
for i=1:Ny
    ypositions=repmat('r',1,i);
    Etensor_r(:,1,1,i)=h5read(strcat('exxkx',num2str(kx),ypositions,'.h5'),'//ex').*w;
    %Etensor_i(:,1,1,i)=h5read(strcat('exxkx',num2str(kx),ypositions,'.h5'),'//ex.i').*w;
    Etensor_r(:,1,2,i)=h5read(strcat('exykx',num2str(kx),ypositions,'.h5'),'//ex').*w;
    %Etensor_i(:,1,2,i)=h5read(strcat('exykx',num2str(kx),ypositions,'.h5'),'//ex.i').*w;
    Etensor_r(:,1,3,i)=h5read(strcat('exzkx',num2str(kx),ypositions,'.h5'),'//ex').*w;
    %Etensor_i(:,1,3,i)=h5read(strcat('exzkx',num2str(kx),ypositions,'.h5'),'//ex.i').*w;
    Etensor_r(:,2,1,i)=h5read(strcat('eyxkx',num2str(kx),ypositions,'.h5'),'//ey').*w;
    %Etensor_i(:,2,1,i)=h5read(strcat('eyxkx',num2str(kx),ypositions,'.h5'),'//ey.i').*w;
    Etensor_r(:,2,2,i)=h5read(strcat('eyykx',num2str(kx),ypositions,'.h5'),'//ey').*w;
    %Etensor_i(:,2,2,i)=h5read(strcat('eyykx',num2str(kx),ypositions,'.h5'),'//ey.i').*w;
    Etensor_r(:,2,3,i)=h5read(strcat('eyzkx',num2str(kx),ypositions,'.h5'),'//ey').*w;
    %Etensor_i(:,2,3,i)=h5read(strcat('eyzkx',num2str(kx),ypositions,'.h5'),'//ey.i').*w;
    Etensor_r(:,3,1,i)=h5read(strcat('ezxkx',num2str(kx),ypositions,'.h5'),'//ez').*w;
    %Etensor_i(:,3,1,i)=h5read(strcat('ezxkx',num2str(kx),ypositions,'.h5'),'//ez.i').*w;
    Etensor_r(:,3,2,i)=h5read(strcat('ezykx',num2str(kx),ypositions,'.h5'),'//ez').*w;
    %Etensor_i(:,3,2,i)=h5read(strcat('ezykx',num2str(kx),ypositions,'.h5'),'//ez.i').*w;
    Etensor_r(:,3,3,i)=h5read(strcat('ezzkx',num2str(kx),ypositions,'.h5'),'//ez').*w;
    %Etensor_i(:,3,3,i)=h5read(strcat('ezzkx',num2str(kx),ypositions,'.h5'),'//ez.i').*w;
end
Etensor=Etensor_r;%+1i*Etensor_i;
clear Etensor_r Etensor_i 
