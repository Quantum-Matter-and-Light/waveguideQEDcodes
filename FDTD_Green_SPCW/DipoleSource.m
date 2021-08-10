%%
function DipoleSource(h,cont,x,y,z,freq,bw,thetad,phid)
lcommand(h,'adddipole',cont);
%lset(h,cont,'name',['DSource at (x,y,z)=' num2str([x y z]) ' with freq=' num2str(freq) ' bandwidth=' num2str(bw)]);
lset(h,cont,'name','DSource');
lset(h,cont,'x',x);lset(h,cont,'y',y);lset(h,cont,'z',z);
lset(h,cont,'set frequency',1);
lset(h,cont,'override global source settings',true);
lset(h,cont,'center frequency',freq);lset(h,cont,'frequency span',bw);
lset(h,cont,'optimize for short pulse',true);
lset(h,cont,'theta',thetad);
lset(h,cont,'phi',phid);

end