function  setglobal(name, val)
%SETGLOBAL Summary of this function goes here
%   Detailed explanation goes here

eval(['global ' name ';']);
eval([name '= val;']);

end

