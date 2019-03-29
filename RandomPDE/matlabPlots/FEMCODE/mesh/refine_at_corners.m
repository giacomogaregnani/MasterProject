function mesh = refine_at_corners(mesh, N)
for i = 1:N
    bar = get_rc(mesh);
    d1 = bsxfun(@minus, bar, [0,0]); d1 = sqrt(d1(:,1).^2 + d1(:,2).^2);
    d2 = bsxfun(@minus, bar, [0,1]); d2 = sqrt(d2(:,1).^2 + d2(:,2).^2);
    d3 = bsxfun(@minus, bar, [1,1]); d3 = sqrt(d3(:,1).^2 + d3(:,2).^2);
    d4 = bsxfun(@minus, bar, [1,0]); d4 = sqrt(d4(:,1).^2 + d4(:,2).^2);
    [m1,n1] = sort(d1);
    [m2,n2] = sort(d2);
    [m3,n3] = sort(d3);
    [m4,n4] = sort(d4);
    mesh = bisect(mesh,[n1(1:4) n2(1:4) n3(1:4) n4(1:4)]);
end