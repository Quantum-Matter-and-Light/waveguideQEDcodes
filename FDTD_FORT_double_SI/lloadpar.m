function code= lloadpar(h,parname,par,cont)
%path(path, 'C:\Program Files\Lumerical\FDTD\api\matlab')
if isstr(par)
    code = strcat([parname '= ' par ';']);
else
    code=strcat([parname '= ' num2str(par) ';']);
end
disp(code);
if cont 
    appevalscript(h,code);
end
end