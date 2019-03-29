% Plot the grid used for PRE-PROCESSING interpolation

dx = 0.2;
x = -1:dx:1;

figure
plot([-1,-1],[-1,1],'k','LineWidth',2);
hold on
plot([1,1],[-1,1],'k',[-1,1],[-1,-1],'k',[-1,1],[1,1],'k','LineWidth',2)

for i = 2 : length(x) - 1
    plot([-1,1],[x(i),x(i)],'--k')
    plot([x(i),x(i)],[-1,1],'--k')
    plot(x(i)-dx/2 * ones(length(x)-1),x(2)-dx/2:dx:x(end)-dx/2,'or')
end

plot(x(end)-dx/2 * ones(length(x)-1),x(2)-dx/2:dx:x(end)-dx/2,'or')
xlabel('x')
ylabel('y')
axis([-1.1 1.1 -1.1 1.1])
axis equal
