display('permutations and their parity (last column)');
p=perms(1:3);
display([p, permutationparity(p,2)]);