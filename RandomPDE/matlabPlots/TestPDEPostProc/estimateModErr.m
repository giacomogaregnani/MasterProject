function [mE, vE] = estimateModErr(meshRef, mesh, data, xObs, mPrior, sPrior, nModError)

nObs = length(xObs);
e = zeros(nObs, nModError);
for j = 1 : nModError
    t = mPrior + sPrior * randn(size(mPrior));
    kModError = buildField_Andrea(t);
    uRefModErr = solveFwdProblem_Cont(meshRef, kModError, data.f, data.rBC);
    uRefModErr = pointEval(uRefModErr, xObs, meshRef);
    uModErr = solveFwdProblem_Cont(mesh, kModError, data.f, data.rBC);
    uModErr = pointEval(uModErr, xObs, mesh);
    e(:, j) =  uRefModErr - uModErr;
end
mE = mean(e, 2);
vE = cov(e');