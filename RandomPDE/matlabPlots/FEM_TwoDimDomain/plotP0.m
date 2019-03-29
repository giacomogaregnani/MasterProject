function plotP0(mesh, u)
% Plot a piece-wise constant solution u

figure
p = mesh.vertices;
t = mesh.elements(1:3, :);
x = p(1,:);
y = p(2,:);
P = [x(t(:));y(t(:))];
T = reshape(1:size(P,2),[3 size(P,2)/3]);
tmp = [u';u';u'];
trisurf(T',P(1,:),P(2,:),tmp(:))
colorbar
view(2)