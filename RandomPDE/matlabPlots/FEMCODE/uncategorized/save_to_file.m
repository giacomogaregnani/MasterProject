function save_to_file(file, varargin)
%SAVE_TO_FILE Summary of this function goes here
%   Detailed explanation goes here
N = numel(varargin)/2;
for i=1:N
    eval([varargin{2*i-1} '= varargin{' num2str(2*i) '};']);
end
str= ['save(''' file ''',' ];
for i=1:N
    str = [str '''' varargin{2*i-1}  ''''];
    if i<N
        str = [str ','];
    else
        str = [str ');'];
    end
end
eval(str);
end

