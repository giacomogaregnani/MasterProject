clear; clc; close all
load('order2')
errtwo = err;
load('order3')
errthree = err;
load('order4')
errfour = err;
load('order4two')
errfourtwo = err;

loglog(h, errtwo, '-o')
hold on
loglog(h, errthree, '+-')
loglog(h, errfour, '<-')
loglog(h, errfourtwo, '*-') 
loglog(h, 10 * h.^2, 'k--')
loglog(h, 10 * h.^3, 'k-.')
loglog(h, 10 * h.^4, 'k:')

xlabel('$h$')
ylabel('error')
legend('$p = 1$', '$p = 1.5$', '$p = 2$', '$p = 2.5$', 'slope 2', ...
    'slope 3', 'slope 4', 'Location', 'SE')  
