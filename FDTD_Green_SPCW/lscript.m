function code=lscript(varargin)
code='';
for i1=1:nargin
    code=[code, ...
        varargin{i1} ';'];
end
end