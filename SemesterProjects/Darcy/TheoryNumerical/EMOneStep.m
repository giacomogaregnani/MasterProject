function X = EMOneStep(X,f,g,dW,h)
% One STEP of EM 

X = X + f(X(1),X(2))*h + g(X(1),X(2))*dW;

end