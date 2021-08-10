function [g]=makegt(Courant,res,endtime,trans_freq,Sigma)
%trans_freq=1.174;
%Courant=0.25;
%res=50;
%endtime=10;
    Resolution=res;
    dt=Courant/Resolution;
    Tfft=100;
    if Tfft < endtime
        Tfft=endtime*100;
    elseif Tfft < 1000
        Tfft=1000;
    end
    t=0:Courant/Resolution:endtime;
    Nfft=round(Tfft/dt);
    dg=zeros(Nfft,1);
    for ii=1:Nfft/2
        xi=2*pi*ii/(Nfft*dt);
        dg(ii)=-1i*xi*sqrt(1+1i*Sigma/xi)*(1+(1i*0.5*Sigma/xi))/(xi*sqrt(1+1i*Sigma/xi)+(trans_freq));
    end    

    for ll=1:Nfft/2
        xi=2*pi*ll/(Nfft*dt);
        ddg(ll)=dg(ll)-(-1i+0.5.*sqrt((1i.*Sigma.*Sigma.*Sigma)/xi)/(trans_freq));
    end    
    
    ggt=zeros(Nfft,1);
    ggt=ifft(ddg,Nfft);
    gtt=ifft(dg,Nfft);
    Nt=round(endtime/dt);
    dxi=1/(dt);                     %Note that FFTW normalization is different from the one in Matlab
    for jj=1:Nt
        t=jj*dt;
        g(jj)=2*((ggt(jj)*dxi)+1i/(2*pi*t)+1i*0.25*sqrt(Sigma*Sigma*Sigma/(t*pi))/(trans_freq));
    end
    for jj=1:Nt
        t=jj*dt;
        ggg(jj)=2*(gtt(jj)*dxi);
    end
    
