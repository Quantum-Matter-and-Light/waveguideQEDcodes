function [data,bloc]=readBands(file,savename)
if ~nargin | isempty(file) | ~ischar(file)
    error('Unknown or undefined argument')
end
if ~nargout
    error('Unknown return value.');
end;
[fid,msg]=fopen(file,'rt');
if fid<0
    error(msg);
end;
data={};
line=fgetl(fid);
    bloc.file=fopen(fid);
    bloc.grid=zeros(3,1);
    bloc.size=zeros(3,1);
    bloc.lattice=zeros(3);
    bloc.reciprocal=zeros(3);
    bloc.geometry={};
    bloc.polarity='';
    bloc.vectors=[];
    bloc.bands=[];
    bloc.ranges=[];
    bloc.xvelocity=[];
    bloc.yvelocity=[];
    bloc.zvelocity=[];
    bloc.yparity=[];
    bloc.zparity=[];
    bloc.time=0;
while ischar(line)

    while ~length(strmatch('Grid size',line))
        line=fgetl(fid);
        if ~ischar(line)
            return;
        end;
    end;
    bloc.grid=sscanf(line,'Grid size is %d x %d x %d');
    line=fgetl(fid);
    while ischar(line) & ~length(strmatch('Lattice vectors',line))
        line=fgetl(fid);
    end;
    for i=1:3
        line=fgetl(fid);
        if ~ischar(line)
            'Corrupt lattice definition.'
            return;
        end;
        bloc.lattice(i,:)=sscanf(line(findstr(line,'('):length(line)),'(%f, %f, %f)')';
        bloc.size(i)=norm(bloc.lattice(i,:));
    end;
    line=fgetl(fid);
    while ischar(line) & ~length(strmatch('Reciprocal lattice vectors',line))
        line=fgetl(fid);
    end;
    for i=1:3
        line=fgetl(fid);
        if ~ischar(line)
            'Corrupt reciprocal lattice definition.'
            return;
        end;
        bloc.reciprocal(i,:)=sscanf(line(findstr(line,'('):length(line)),'(%f, %f, %f)')';
    end;
    line=fgetl(fid);
    while ischar(line) & ~length(strmatch('Geometric objects',line))
        line=fgetl(fid);
    end;
    line=fgetl(fid);
    while ischar(line) & ~length(strmatch('Geometric object t',line))
        bloc.geometry{length(bloc.geometry)+1}=line;
        line=fgetl(fid);
    end;
    i=0;
    n=[];
    line=fgetl(fid);
    while ischar(line) & isempty(n)
        n=sscanf(line,'%d k-points:');
        line=fgetl(fid);
    end;
    while ischar(line) & ~length(strmatch('Solving for band p',line))
        line=fgetl(fid);
    end;
    if ~ischar(line)
        'Corrupt polarity definition.'
        return;
    end;
    bloc.polarity=line(findstr(line,': ')+2:length(line)-1);
    line=fgetl(fid);
    while ischar(line) & ~length(strmatch('Band ',line))
        s=findstr(line,':, ');
        if length(s) ==1 & s>5
            switch line(s-9:s-1)
                case 'evenfreqs'        % band frequency
                    if i
                        v=sscanf(line(s+1:length(line)),', %f',inf)';
                        if v(1) ~= i
                            'Synchronization lost.'
                        end;
                        bloc.bands(i,:)=v(6:length(v));
                        bloc.vectors(i,:)=v(2:5);
                    else
                        bloc.bands=zeros(n,length(findstr(line,','))-5);
                    end;
                    i=i+1
                case 'nvelocity'        % x/y/z group velocity
                    s=findstr(line,',');
                    if sscanf(line(s(1):s(2)),', %d') ~= i-1
                        'Synchronization lost.'
                    end;
                    v=sscanf(line(s(2):length(line)),', #(%f %f %f)',[3 inf]);
                    if i==2
                        bloc.xvelocity=zeros(n,size(v,2));
                        bloc.yvelocity=bloc.xvelocity;
                        bloc.zvelocity=bloc.xvelocity;
                    end;
                    bloc.xvelocity(i-1,:)=v(1,:);
                    bloc.yvelocity(i-1,:)=v(2,:);
                    bloc.zvelocity(i-1,:)=v(3,:);
                    
            end;
        end;
        line=fgetl(fid);
    end;
    n=size(bloc.bands,2);     % band limits
    bloc.ranges=zeros(2,n);
    for i=1:n
        if ~ischar(line)
            return;
        end;
        if sscanf(line(6:10),'%d') ~= i
            'Synchronization lost.'
        end;
        bloc.ranges(:,i)=sscanf(line(findstr(line,': ')+1:length(line)),'%f at #(%*f %*f %*f) to %f')';
        line=fgetl(fid);
    end;
    while ischar(line)           % total elapsed time
        if  strcmp(line(1:5),'total')
            line=line(findstr(line,':')+2:length(line));
            n=findstr(line,',');
            while length(n)
                bloc.time=bloc.time*60+sscanf(line,'%d');
                line=line(n(1)+2:length(line));
                n=findstr(line,',');
            end;
            bloc.time=bloc.time*60+sscanf(line,'%d');
            break;
        end;
        line=fgetl(fid);
    end;
    data{length(data)+1}=bloc;
    line=fgetl(fid);
end;
fclose(fid);
bloc={};
