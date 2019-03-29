function varargout = load_from_file(file, varargin)
%LOAD_FROM_FILE Summary of this function goes here
%   Detailed explanation goes here
N = numel(varargin);
str= ['load(''' file ''',' ];
for i=1:N
    str = [str '''' varargin{i}  ''''];
    if i<N
        str = [str ','];
    else
        str = [str ');'];
    end
end

eval(str);
for i=1:N
    eval(['varargout{' num2str(i) '} = ' varargin{i} ';']);
end
end

