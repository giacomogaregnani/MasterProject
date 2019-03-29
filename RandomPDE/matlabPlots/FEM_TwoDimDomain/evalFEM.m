function uEval = evalFEM(MESH, U, X)
% X has to be a matrix with 2 rows and n columns where n is the number of
% points where the solution has to be evaluated

uInt = pdeInterpolant(MESH.nodes, MESH.elements, U);
uEval = evaluate(uInt, X);

end