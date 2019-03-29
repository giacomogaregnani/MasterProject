function X = EMOneStep(X,f,g,dW,h)
% One step of EULER MARUYAMA

X = X + f(X)*h + g(X)*dW;

end