function flag = get_bdflag(mesh, fun)
%GET_BDFLAG Summary of this function goes here
%   Detailed explanation goes here

[subsim,isbdsubsim,subsim2elem] = ...
	auxstructure(mesh,'subsim','isbdsubsim','subsim2elem');
flagged = fun(get_rc(mesh,subsim(isbdsubsim,:)));
isbdsubsim(isbdsubsim ) = flagged;
flag = subsim2elem(isbdsubsim,[1 3]);
end

