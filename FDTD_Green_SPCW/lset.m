function code=lset(h,cont,varargin)
code='';
%nargin
for i1=1:nargin-2
    if isstr(varargin{i1})
        code=[code '"' varargin{i1} '"'];
    else
        code=[code '' num2str(varargin{i1}) ''];
    end
    if i1==nargin-2
        code=[code];
    else
        code=[code ','];
    end
end
code=['set(' code ')'];
lcommand(h,code,cont);
end