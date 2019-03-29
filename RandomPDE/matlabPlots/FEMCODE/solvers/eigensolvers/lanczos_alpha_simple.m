function [eigmu, eiglbmu, eigvecarray, niter] = ...
    lanczos_alpha_simple(A, Chinv, ChinvT, n)
% We solve for the first eigenvalue of  
%   [bch'\amat*bch\]*eigvec = eigval*eigvec 
% which is equivalent to
%   amat*eigvec = eigval*[bch'*bch]*eigvec 
% The eigenvector (eigvec) is then given in the second format

nitermx = 500;
tol = 1e-10;

binit = rand(n,1);  % initial vector
opts.disp = 0;

%% LANCZOS ALGORITHM
% q is the set of (orthogonal) Lanczos vectors./
q = zeros(n,nitermx);
tri = sparse(nitermx,nitermx);
betaim1 = 0;
qim1 = zeros(n,1);
normb = sqrt(binit'*binit);
qi = binit/normb;
q(:,1) = qi;
niter = nitermx;
displaystr = '';

for i = 1:nitermx;
    % Find next Lanczos vector and next row of Lanczos tridiagonal matrix:
    v = ChinvT(A(Chinv(qi))) - betaim1*qim1;
    alphai = qi'*v;
    v = v  - alphai*qi;
    betai = sqrt(v'*v);
    tri(i,i) = alphai;
    if(i ~= 1);
        tri(i,i-1) = betaim1;
    end
    if(i ~= nitermx);
        tri(i,i+1) = betai;
    end
    betaim1 = betai;
    qim1 = qi;
    qi = v/betai;
    if(i ~= nitermx)
        q(:,i+1) = qi;
    end
    if i <= 1
        continue;
    end
    % Compute eigenvalues of current (i x i) Lanczos tridiagonal "stiffness" matrx tri:
    [eigvecmu,eigmu] = eigs(tri(1:i,1:i), eye(i), 1, 'SA', opts);
    % Compute error bound errbnd for Ritz values relative to nearest exact
    % (truth) eigenvalue (see Demmel)
    errbnd = betai*abs(eigvecmu(i));
    eiglbmu = eigmu  - errbnd;
    ratio = errbnd/(tol*eigmu/(1 + tol));
    eigvec = q(:,1:size(eigvecmu,1))*eigvecmu;     % Get the true eigvector
    
    % PRINT OUTPUTe
    if ~strcmp(displaystr,'')
        fprintf(repmat('\b',[1,numel(displaystr)]));
    end
    displaystr = sprintf('LANCZOS:: size %d, iter: %d, eigval: %f, ratio: %f',n,i,eigmu,ratio);   
    fprintf(displaystr);   
    
    % Terminate if ratio less than unity for all eigenvalues:
    if abs(ratio) <= 1;
        niter = i;
        break;
    end   
end

fprintf('\n');

%% OUTPUT MANIPULATION
% be careful as the two eigenproblems are not the same!
eigvecarray = Chinv(eigvec(:,1));
end