function X = EMOneStep(X,f,g,dW,h)
% One step of EULER MARUYAMA to integrate with DEM and CEM

X = X + f(X(1),X(2))*h + g(X(1),X(2))*dW;

end