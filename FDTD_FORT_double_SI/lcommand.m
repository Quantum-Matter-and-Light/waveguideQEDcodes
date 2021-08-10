function lcommand(h,str,cont)
code=strcat([str ';']);
disp(code);

if cont
    appevalscript(h,code);
end
end