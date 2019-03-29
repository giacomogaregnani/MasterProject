function X = EMOneStepDarcy(X,u,sigma,dW,h)
% One STEP of EM in the DARCY case

X = X + u*h + sigma*dW;

end