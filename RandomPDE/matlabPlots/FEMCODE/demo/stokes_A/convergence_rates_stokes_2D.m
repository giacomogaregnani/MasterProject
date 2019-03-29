clear vp sol;
vp.a = 1;
vp.f = @fstokesex;
uelemtype = {'p1b', 'p2', 'p3'};
pelemtype = {'p1', 'p1', 'p2'};
E = numel(uelemtype);
vp.bc = 'zero_dirichlet';


Nmax = 4;
[h1,l2,h1s,dof] = deal(zeros(Nmax,3));
[usol, ufemspace, psol, pfemspace, mesh, vps] = deal(cell(E,Nmax));

for j=1:E
    vp.elemtype = uelemtype{j};
    vp.pelemtype = pelemtype{j};
    for i=1:Nmax
        thismesh = structured_mesh([0,1,0,1], 2^i, struct('centre',true));
        [usol{i,j}, ufemspace{i,j}, psol{i,j}, pfemspace{i,j}, ...
            mesh{i,j}, vps{i,j}] = stokes(thismesh, vp);
        dof(i,j) = 2*ufemspace{i,j}.ndof + pfemspace{i,j}.ndof;
        [h1(i,j), l2(i,j), h1s(i,j)] = ...
            get_H1error_exact(mesh{i,j}, ufemspace{i,j}, usol{i,j}, @ex_stokes);
    end
end

figure;
loglog(dof(:,1),l2(:,1),'-x', ...
    dof(:,2),l2(:,2),'-x', ...
    dof(:,3),l2(:,3),'-x');
legend('p1b','p2','p3');


