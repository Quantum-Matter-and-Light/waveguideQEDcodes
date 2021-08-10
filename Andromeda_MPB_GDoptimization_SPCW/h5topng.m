function []=h5topng(fname,rptnum,kk,nband,res)

Eps_zz=h5read(strcat(fname,'-epsilon.h5'),'//epsilon.zz');
if kk<10
    if nband<10
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-epsilon.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.x.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.y.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.z.zeven.h5'])
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -z ' num2str(size(Eps_zz,1)/2) ' -d data-new ' fname '-epsilon.h5'])
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqR0' num2str(nband) '.h5 ' fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.x.zeven.h5:x.r-new ' ...
            fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.y.zeven.h5:y.r-new ' fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.z.zeven.h5:z.r-new']);
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqI0' num2str(nband) '.h5 ' fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.x.zeven.h5:x.i-new ' ...
            fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.y.zeven.h5:y.i-new ' fname '-e.k0' num2str(kk) '.b0' num2str(nband) '.z.zeven.h5:z.i-new']);
        unix(['h5math -e "d1+d2" eSqb0' num2str(nband) '.h5 eSqR0' num2str(nband) '.h5 eSqI0' num2str(nband) '.h5']);
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -x ' num2str(size(Eps_zz,3)/2+res) ...
            ' -o eSqb0' num2str(nband) '0' num2str(kk) '.xslice.png eSqb0' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -y ' num2str(size(Eps_zz,2)/2) ...
            ' -o eSqb0' num2str(nband) '0' num2str(kk) '.yslice.png eSqb0' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -z ' num2str(size(Eps_zz,1)/2) ...
            ' -o eSqb0' num2str(nband) '0' num2str(kk) '.zslice.png eSqb0' num2str(nband) '.h5']);
    elseif nband>10
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-epsilon.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k0' num2str(kk) '.b' num2str(nband) '.x.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k0' num2str(kk) '.b' num2str(nband) '.y.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k0' num2str(kk) '.b' num2str(nband) '.z.zeven.h5'])
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -z ' num2str(size(Eps_zz,1)/2) ' -d data-new ' fname '-epsilon.h5'])
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqR' num2str(nband) '.h5 ' fname '-e.k0' num2str(kk) '.b' num2str(nband) '.x.zeven.h5:x.r-new ' ...
            fname '-e.k0' num2str(kk) '.b' num2str(nband) '.y.zeven.h5:y.r-new ' fname '-e.k0' num2str(kk) '.b' num2str(nband) '.z.zeven.h5:z.r-new']);
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqI' num2str(nband) '.h5 ' fname '-e.k0' num2str(kk) '.b' num2str(nband) '.x.zeven.h5:x.i-new ' ...
            fname '-e.k0' num2str(kk) '.b' num2str(nband) '.y.zeven.h5:y.i-new ' fname '-e.k0' num2str(kk) '.b' num2str(nband) '.z.zeven.h5:z.i-new']);
        unix(['h5math -e "d1+d2" eSqb' num2str(nband) '.h5 eSqR' num2str(nband) '.h5 eSqI' num2str(nband) '.h5']);
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -x ' num2str(size(Eps_zz,3)/2+res) ...
            ' -o eSqb' num2str(nband) '0' num2str(kk) '.xslice.png eSqb' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -y ' num2str(size(Eps_zz,2)/2) ...
            ' -o eSqb' num2str(nband) '0' num2str(kk) '.yslice.png eSqb' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -z ' num2str(size(Eps_zz,1)/2) ...
            ' -o eSqb' num2str(nband) '0' num2str(kk) '.zslice.png eSqb' num2str(nband) '.h5']);
    end
elseif kk>10
    if nband<10
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-epsilon.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k' num2str(kk) '.b0' num2str(nband) '.x.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k' num2str(kk) '.b0' num2str(nband) '.y.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k' num2str(kk) '.b0' num2str(nband) '.z.zeven.h5'])
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -z ' num2str(size(Eps_zz,1)/2) ' -d data-new ' fname '-epsilon.h5'])
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqR0' num2str(nband) '.h5 ' fname '-e.k' num2str(kk) '.b0' num2str(nband) '.x.zeven.h5:x.r-new ' ...
            fname '-e.k' num2str(kk) '.b0' num2str(nband) '.y.zeven.h5:y.r-new ' fname '-e.k' num2str(kk) '.b0' num2str(nband) '.z.zeven.h5:z.r-new']);
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqI0' num2str(nband) '.h5 ' fname '-e.k' num2str(kk) '.b0' num2str(nband) '.x.zeven.h5:x.i-new ' ...
            fname '-e.k' num2str(kk) '.b0' num2str(nband) '.y.zeven.h5:y.i-new ' fname '-e.k' num2str(kk) '.b0' num2str(nband) '.z.zeven.h5:z.i-new']);
        unix(['h5math -e "d1+d2" eSqb0' num2str(nband) '.h5 eSqR0' num2str(nband) '.h5 eSqI0' num2str(nband) '.h5']);
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -x ' num2str(size(Eps_zz,3)/2+res) ...
            ' -o eSqb0' num2str(nband) num2str(kk) '.xslice.png eSqb0' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -y ' num2str(size(Eps_zz,2)/2) ...
            ' -o eSqb0' num2str(nband) num2str(kk) '.yslice.png eSqb0' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -z ' num2str(size(Eps_zz,1)/2) ...
            ' -o eSqb0' num2str(nband) num2str(kk) '.zslice.png eSqb0' num2str(nband) '.h5']);
        
    elseif nband>10
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-epsilon.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k' num2str(kk) '.b' num2str(nband) '.x.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k' num2str(kk) '.b' num2str(nband) '.y.zeven.h5'])
        unix(['~/apps/mpb-1.5/bin/mpb-data -y ' num2str(rptnum) ' -r -n ' num2str(res) ' ' fname '-e.k' num2str(kk) '.b' num2str(nband) '.z.zeven.h5'])
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -z ' num2str(size(Eps_zz,1)/2) ' -d data-new ' fname '-epsilon.h5'])
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqR' num2str(nband) '.h5 ' fname '-e.k' num2str(kk) '.b' num2str(nband) '.x.zeven.h5:x.r-new ' ...
            fname '-e.k' num2str(kk) '.b' num2str(nband) '.y.zeven.h5:y.r-new ' fname '-e.k' num2str(kk) '.b' num2str(nband) '.z.zeven.h5:z.r-new']);
        unix(['h5math -e "d1^2+d2^2+d3^2" eSqI' num2str(nband) '.h5 ' fname '-e.k' num2str(kk) '.b' num2str(nband) '.x.zeven.h5:x.i-new ' ...
            fname '-e.k' num2str(kk) '.b' num2str(nband) '.y.zeven.h5:y.i-new ' fname '-e.k' num2str(kk) '.b' num2str(nband) '.z.zeven.h5:z.i-new']);
        unix(['h5math -e "d1+d2" eSqb' num2str(nband) '.h5 eSqR' num2str(nband) '.h5 eSqI' num2str(nband) '.h5']);
        
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -x ' num2str(size(Eps_zz,3)/2+res) ...
            ' -o eSqb' num2str(nband) num2str(kk) '.xslice.png eSqb' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -y ' num2str(size(Eps_zz,2)/2) ...
            ' -o eSqb' num2str(nband) num2str(kk) '.yslice.png eSqb' num2str(nband) '.h5']);
        unix(['~/apps/h5utils-1.13/bin/h5topng -C ' fname '-epsilon.h5:data-new -c bluered -Z -z ' num2str(size(Eps_zz,1)/2) ...
            ' -o eSqb' num2str(nband) num2str(kk) '.zslice.png eSqb' num2str(nband) '.h5']);
        
    end
end

