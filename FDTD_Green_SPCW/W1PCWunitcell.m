function W1PCWunitcell(h,cont,nSiN,tSiN,x0,y0,z0,a0,wid,thk,ny,r1,r2,r3,s1,s2,offset,Alpha);
lcommand(h,'addstructuregroup',cont);
lset(h,cont,'name','W1PhCW_Structure');
lcommand(h,'addstructuregroup',cont);
lset(h,cont,'name','W1PhCW_cell');

wX=a0;
L1=s1+offset+(sqrt(3)/2)*a0;
L2=s2+offset+sqrt(3)*a0;
wY=a0*sqrt(3)*ny+2*L2;
elpa=r1*Alpha;
elpb=r1/Alpha;

lcommand(h,'addrect',cont);
lset(h,cont,'name',['W1PCW slab' num2str(round(x0/a0))]);
lset(h,cont,'index',nSiN);
lset(h,cont,'x',x0);
lset(h,cont,'y',y0);
lset(h,cont,'z',z0);
lset(h,cont,'x span',wX);
lset(h,cont,'y span',wY);
lset(h,cont,'z span',thk);

lcommand(h,'addrect',cont);
lset(h,cont,'name','W1PCW slot');
lset(h,cont,'index',1);
lset(h,cont,'x',x0);
lset(h,cont,'y',y0);
lset(h,cont,'z',z0);
lset(h,cont,'x span',a0);
lset(h,cont,'y span',wid);
lset(h,cont,'z span',thk);
sign=[-1 1];

for i=1:2;
    for j=1:2;
        lcommand(h,'addcustom',cont);
        lset(h,cont,'name','W1PCW squircle');
        lset(h,cont,'create 3D object by','extrusion');
        lset(h,cont,'equation units','m');
        lset(h,cont,'index',1);
        lset(h,cont,'make nonsymmetric',1);
        lset(h,cont,'equation 1',[num2str(elpb) '*(' num2str(1) '-(x/' num2str(elpa) ')^' num2str(4) ')^' num2str(0.25)]);
        lset(h,cont,'equation 2',[num2str(elpb) '*(' num2str(1) '-(x/' num2str(elpa) ')^' num2str(4) ')^' num2str(0.25)]);
        lset(h,cont,'x',x0+sign(i)*a0/2);
        lset(h,cont,'x span',2*elpa);
        lset(h,cont,'y',y0+sign(j)*L1);
        lset(h,cont,'y span',2*elpb); % computational size
        lset(h,cont,'z',z0);
        lset(h,cont,'z span',thk);
        lset(h,cont,'detail',1)
    end
end
lcommand(h,'addcircle',cont);
lset(h,cont,'name','W1PCW 2nd hole up');
lset(h,cont,'index',1);
lset(h,cont,'x',x0);
lset(h,cont,'y',y0+L2);
lset(h,cont,'z',z0);
lset(h,cont,'z span',thk);
lset(h,cont,'radius',r2);

lcommand(h,'addcircle',cont)
lset(h,cont,'name','W1PCW 2nd hole down');
lset(h,cont,'index',1);
lset(h,cont,'x',x0);
lset(h,cont,'y',y0-L2);
lset(h,cont,'z',z0);
lset(h,cont,'z span',thk);
lset(h,cont,'radius',r2);

for k=1:2
    for j=1:2
        for i=1:ny/2
            lcommand(h,'addcircle',cont)
            lset(h,cont,'name',['W1PCW other holes1' num2str(i) num2str(j) num2str(k)]);
            lset(h,cont,'index',1);
            lset(h,cont,'x',x0+(sign(j)*a0/2));
            lset(h,cont,'y',y0+sign(k)*(L2+(-sqrt(3)/2)*a0+i*sqrt(3)*a0));
            lset(h,cont,'z',z0);
            lset(h,cont,'z span',thk);
            lset(h,cont,'radius',r3);
        end
    end
end
for k=1:2
    for i=1:ny/2
        lcommand(h,'addcircle',cont)
        lset(h,cont,'name',['W1PCW other holes2' num2str(i) num2str(j) num2str(k)]);
        lset(h,cont,'index',1);
        lset(h,cont,'x',x0);
        lset(h,cont,'y',y0+sign(k)*(L2+(i*sqrt(3)*a0)));
        lset(h,cont,'z',z0);
        lset(h,cont,'z span',thk);
        lset(h,cont,'radius',r3);
    end
end
N1 = 30; %% changed to get 50 unit cell // Kyung
ddx=x0-a0/2;
object_counter=1;
lcommand(h,'selectpartial("W1PCW")',true);
lcommand(h,'addtogroup("W1PhCW_cell")',true);
lcommand(h,'select("W1PhCW_cell")',true);
for i=1:N1 
    if (object_counter > 1) 
       x1 = 1*a0*(i-1)+a0/2+ddx;
       lcommand(h,strcat('copy(',num2str(x1-x0),',',num2str(0),',',num2str(0),')'),true);
       x0 = x1;
    end
    object_counter = object_counter + 1;
end
lcommand(h,'select("W1PhCW_cell")',true);
lcommand(h,'addtogroup("W1PhCW_Structure")',true);
lcommand(h,'select("W1PhCW_Structure")',true);
lcommand(h,'redrawoff',true);
object_counter=1;
x0=ddx;
N2=2;
for i=1:N2 
    if (object_counter > 1) 
       x1 = N1*a0*(i-1)+ddx;
       lcommand(h,strcat('copy(',num2str(-x1+x0),',',num2str(0),',',num2str(0),')'),true);
       x0 = x1;
    end
    object_counter = object_counter + 1;
end
lcommand(h,'unselectall',true);
lcommand(h,'redrawon',true);

end