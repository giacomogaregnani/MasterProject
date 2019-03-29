function l = likelihoodFEM(param, field, f, x, xObs, observations, sigma)

    F = assembleRHS(f, x);
    A = assembleMatrix(@(xx) field(xx, param), x);
    uVec = [0; A \ F; 0];
  
    u = interp1(x, uVec, xObs);
    
    diff = observations - u;
    l = -1 / (2 * sigma^2) * (diff') * diff;
     
end