function expr = evalmu(expr, mu)
evalstr = strrep(strrep(strrep(char(expr),'*','.*'),'^','.^'),'/','./');
for ii=1:size(mu,2)
    evalstr = strrep(evalstr, ...
        ['mu' num2str(ii)], ...
        ['mu(:,' num2str(ii) ')']);
end
eval(['expr = ' evalstr ';']);
end