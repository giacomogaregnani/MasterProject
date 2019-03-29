%% PERFORMANCE TEST

node2 = [0,0; 0,2; 1,2; 1,1; 2,1; 2,0];
bd2 = [1,2; 2,3; 3,4; 4,5; 5,6; 6,1];
node3 = [0,0,0; 1,0,0; 0,1,0; -1,-1,1; -1,-1,-1];
bd3 = [1,2,4; 2,3,4; 3,1,4; 2,1,5; 3,2,5; 1,3,5];

display('500 000 points in 2D - Lshape');
N=500000;
p2=rand(N,2)*3-0.5;
tic;
dpoly( p2, node2, bd2, options );
toc;


display('500 000 points in 3D - simple nonconvex');
p3=rand(N,3)*3-1.5;
tic;
dpoly( p3, node3, bd3, options );
toc;

options.inside = 0;


display('500 000 points in 2D - Lshape, no inside');
p2=rand(N,2)*3-0.5;
tic;
dpoly( p2, node2, bd2, options );
toc;

display('500 000 points in 3D - simple nonconvex, no inside');
p3=rand(N,3)*3-1.5;
tic;
dpoly( p3, node3, bd3, options );
toc;