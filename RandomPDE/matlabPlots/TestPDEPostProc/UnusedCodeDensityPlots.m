if plotDensities
        if exSol1p
            
            % Build exact distribution
            figure
            hold on

            uEx = @(theta, x) exp(-theta) * (2 * pi)^(-2) * sin(2 * pi * x);
%             uEx = @(theta, x) exp(-theta) * 2 / 3 * (x - x.^3);
            
            post = @(theta) exp(-0.5 * theta.^2 - 0.5 / obsNoise^2 * (uObs - uEx(theta, xObs)).^2);
            logPost = @(theta) -0.5 * theta.^2 - 0.5 / obsNoise^2 * (uObs - uEx(theta, xObs)).^2;
            thetaPlot = linspace(min(thetaAllProbEff) - 0.5, max(thetaAllProbEff) + 0.5, 10000);
            postNorm = post(thetaPlot) / trapz(thetaPlot, post(thetaPlot));
            
            plot(thetaPlot, postNorm);
            % Build posteriors of the chain
            plotDensitiesMCMC(thetaAllEff, thetaRef);
            plotDensitiesMCMC(thetaAllProbEff, thetaRef);
            legend('exact posterior', 'det', 'true', 'prob')
                       
            hPlot = (max(thetaPlot) - min(thetaPlot)) / length(thetaPlot);
            postCDF = cumsum(postNorm) * hPlot;
            [postCDF, mask] = unique(postCDF);
            postCDF(1) = 0;
            postCDF(end) = 1;
            thetaPlot = thetaPlot(mask);            
            sample = interp1(postCDF, thetaPlot, rand(100000, 1));
            [pTex, pUex] = plotMCMCresultsPDE_Andrea(meshRef, sample', thetaRef, f, rBC, uRef, xObs, uObs);
            set(pTex.CurrentAxes, 'yLim', yLim);
            set(pUex.CurrentAxes, 'yLim', yLimU);

            
        elseif exSol2p
            
            thetaSample1 = linspace(thetaRef(1) - 10 * sqrt(obsNoise), thetaRef(1) + 10 * sqrt(obsNoise), 400);
            thetaSample2 = linspace(thetaRef(2) - 10 * sqrt(obsNoise), thetaRef(2) + 10 * sqrt(obsNoise), 400);            
            
            uEx = @(x, k) sin(2 * pi * x) / (k(1) * (2 * pi)^2) .* (x < 0.5) ...
                + sin(2 * pi * x) / (k(2) * (2 * pi)^2) .* (x >= 0.5);
            
            k = 1;
            posterior = zeros(length(thetaSample1));
            for theta1 = thetaSample1
                j = 1;
                for theta2 = thetaSample2
                    theta = [theta1; theta2];
                    uExObs = uEx(xObs, exp(theta))';
                    posterior(k, j) = -1 / (2 * obsNoise^2) * (uExObs - uObs)' * (uExObs - uObs) + ...
                        -1 / 2 * (theta - prior.avg)' * (prior.stddev \ (theta - prior.avg));
                    j = j + 1;
                end
                k = k + 1;
            end
            posterior = exp(posterior);
            posterior = posterior / trapz(thetaSample1, trapz(thetaSample2, posterior));
            
            figure
            hold on
            plotDensitiesMCMC2(thetaAllEff, 'blue');
            plotDensitiesMCMC2(thetaAllProbEff, 'black');
            contour(thetaSample1, thetaSample2, posterior', 10, 'color', 'red')
            plot(thetaRef(1), thetaRef(2), '.r', 'markersize', 20)
            legend('det', 'prob', 'ex', 'true')
            
        else            
            for j = 1 : nParam
                figure
                hold on
                plotDensitiesMCMC(thetaAllEff(j, :), thetaRef(j));
                plotDensitiesMCMC(thetaAllProbEff(j, :), thetaRef(j));
            end
        end
    end
    