function [L,J,F_rec,mf]=spdf(L)
%L=5/2
I=7/2; %for cesium
J1=abs(L-1/2);
J2=abs(L+1/2);
J=J1:1:J2;
F_rec=cell(numel(J),1);
mf=cell(numel(J),1);
for j=1:numel(J);
F=abs(I-J(j)):1:abs(I+J(j));
F_rec{j}=F;
mf_store=cell(numel(F),1);
for k=1:numel(F);
mf_store{k}=-F(k):1:F(k);
end
mf{j}=mf_store;
end

