function [g]=makegt(Courant,res,endtime,trans_freq1,Sigma1,a_meep)
%trans_freq=1.174;
%Courant=0.25;
%res=50;
%endtime=10;
hbar=1.05457148e-34;
c=299792458;
Sigma=Sigma1*c/a_meep;
trans_freq=trans_freq1*c/a_meep;
    Resolution=res;
    dt=(Courant/Resolution);
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
        xi=2*pi*ii/(Nfft*dt*a_meep/c);
        dg(ii)=-1i*xi*sqrt(1+1i*Sigma/xi)*(1+(1i*0.5*Sigma/xi))/(xi*sqrt(1+1i*Sigma/xi)+(trans_freq));
    end
    
    figure(1)
    plot(imag(dg),'b:','linewidth',1.5)     %Test plot: Original g(xi) function
    xlim([0,Nfft/2])
    ylim([-1,1])
    hold on
    grid on
    xlabel('xi')
    ylabel('dg')

    for ll=1:Nfft/2
        xi=2*pi*ll/(Nfft*dt*a_meep/c);
        ddg(ll)=dg(ll)-(-1i+0.5.*sqrt((1i.*Sigma.*Sigma.*Sigma)/xi)/(trans_freq));
    end
    
   plot(imag(ddg),'r-','linewidth',1.5)    %Test plot: after subtraction of divergent terms
   legend('before subtraction','after subtraction')
   title('g(xi/2pi)function')
   hold off
    
    ggt=zeros(Nfft,1);
    ggt=ifft(ddg,Nfft);
    gtt=ifft(dg,Nfft);
    Nt=round(endtime/dt);
    dxi=1/(dt*a_meep/c);           %Note that FFTW normalization is different from the one in Matlab
    for jj=1:Nt
        t=jj*dt*a_meep/c;
        g(jj)=2*((ggt(jj)*dxi)+1i/(2*pi*t)+1i*0.25*sqrt(Sigma*Sigma*Sigma/(t*pi))/(trans_freq));
    end
    for jj=1:Nt
        t=jj*dt*a_meep/c;
        ggg(jj)=2*(gtt(jj)*dxi);
    end
    
    %gt_r=h5read('gxxr.h5','//Gamma.r');
    %gt_i=h5read('gxxr.h5','//Gamma.i');
    
   figure(2)
   plot(imag(ggg),'b:','linewidth',1.5)
   hold on
   plot(imag(g),'r:','linewidth',1.5)
   xlim([0,100])
   legend('without subtract/add procedure','with-')
   title('g(t)function')
   grid on
   xlabel('t')
   ylabel('g')
   hold off