function u = evaluateFEMforError(uInt,x,y,t,param)

u = uInt.evaluate(x,y);
u = reshape(u, length(u)/6, 6);