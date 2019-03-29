function X = EMOneStepDarcy(X,u,sigma,dW,h)
% Performs ONE STEP of EULER MARUYAMA

X = X + u*h + sigma*dW;

end