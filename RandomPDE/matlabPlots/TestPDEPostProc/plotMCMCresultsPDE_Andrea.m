function [pTheta, pU] = plotMCMCresultsPDE_Andrea(meshRef, thetaAll, thetaRef, f, rightBC, uRef, xObs, uObs)

avgKappa = buildField_Andrea(mean(thetaAll, 2));
kappaRef = buildField_Andrea(thetaRef);

pTheta = figure;
hold on
maxLimKappaPrc = buildField_Andrea(prctile(thetaAll', 95));
minLimKappaPrc = buildField_Andrea(prctile(thetaAll', 5));
confInt = [minLimKappaPrc(meshRef.x), fliplr(maxLimKappaPrc(meshRef.x))];
plot(meshRef.x, kappaRef(meshRef.x), 'r');
plot(meshRef.x, avgKappa(meshRef.x), 'b');
fill([meshRef.x, fliplr(meshRef.x)], confInt, 'b', ...
    'linestyle', '--', 'facealpha', 0.1);
legend('true', 'average', 'confidence 0.95', 'Location', 'best')

pU = figure;
close
% avgKappa = buildField_Andrea(mean(thetaAll, 2));
% kappaTheta = buildField_Andrea(thetaRef);
% maxLimKappaPrc = buildField_Andrea(prctile(thetaAll', 95));
% minLimKappaPrc = buildField_Andrea(prctile(thetaAll', 5));
% hold on
% plot(meshRef.x, uRef, 'r')
% uAvg = solveFwdProblem_Cont(meshRef, avgKappa, f, rightBC);
% plot(meshRef.x, uAvg, 'b')
% uMax = solveFwdProblem_Cont(meshRef, maxLimKappaPrc, f, rightBC);
% uMin = solveFwdProblem_Cont(meshRef, minLimKappaPrc, f, rightBC);
% confIntU = [uMin', fliplr(uMax')];
% fill([meshRef.x, fliplr(meshRef.x)], confIntU, 'b', ...
%     'linestyle', '--', 'facealpha', 0.1);
% 
% plot(xObs, uObs, 'or')


end
