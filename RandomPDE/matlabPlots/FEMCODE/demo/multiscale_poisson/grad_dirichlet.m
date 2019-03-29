function out = grad_dirichlet(x, idx)

switch idx
    case 1
        out = [pi*cos(pi*(x(:,1)+x(:,2))); pi*cos(pi*(x(:,1)+x(:,2)))];
    case 2
        out = [-pi*sin(pi*(x(:,1)+x(:,2))); -pi*sin(pi*(x(:,1)+x(:,2)))];
    case 3
        out = 2*[-pi*sin(pi*(x(:,1)+x(:,2))); -pi*sin(pi*(x(:,1)+x(:,2)))];
    case 4
        out = 2*[pi*cos(pi*(x(:,1)+x(:,2))); pi*cos(pi*(x(:,1)+x(:,2)))];
    case 5
        out = 3*[pi*cos(pi*(x(:,1)+x(:,2))); pi*cos(pi*(x(:,1)+x(:,2)))];
    case 6
        out = 3*[-pi*sin(pi*(x(:,1)+x(:,2))); -pi*sin(pi*(x(:,1)+x(:,2)))];
end