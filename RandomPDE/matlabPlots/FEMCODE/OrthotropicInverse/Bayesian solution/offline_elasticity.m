function [AqRB, FqRB, micMesh, micVp, micFemspace, RBspace, RBerror] = offline_elasticity(micMesh, micVp, NRB, NTrain, tolRB, box, Q, parameters_range)

dim = size(micMesh.elem,2) - 1;
SetTrain = 1:NTrain;
XTrain = box(1) + (box(2) - box(1)).*rand(NTrain, dim);
parameters_train = repmat(parameters_range(:,1),1,NTrain) + bsxfun(@times,(parameters_range(:,2) - parameters_range(:,1)),rand(size(parameters_range,1),NTrain));
ChiTrain = cell(3,1);
ChiTrain{1} = XTrain;
ChiTrain{2} = randi(dim+1, NTrain, 1);
ChiTrain{3} = parameters_train;
[Aq, Fq, micMesh, micVp, micFemspace] = get_parameter_independent_forms(micMesh, micVp, Q);
Theta = zeros(Q, NTrain);
for q = 1:Q
    for i = 1:NTrain
        Theta(q,i) = micVp.ThetaField(ChiTrain{1}(i,:), q, ChiTrain{3}(:,i));
    end
end
ND = micFemspace.ndof;
W = get_matrix_scalar_product(micMesh,micVp);
% iW = W\sparse(1:dim*ND, 1:dim*ND, ones(dim*ND,1));

% Compute first basis function
RBspace = zeros(2*ND,1);
while RBspace == zeros(dim*ND,1)
    i1 = randi(NTrain, 1);
    SetTrain = setdiff(SetTrain, i1);
    A = sparse(dim*ND,dim*ND);
    F = sparse(dim*ND, 1);
    for q = 1:Q
        A = A + Theta(q,i1)*Aq{q};
        F = F + Theta(q,i1)*Fq{q}(:,ChiTrain{2}(i1));
    end
    xi1 = A\F;
%     xi1 = reshape(xi1, ND, dim);
%     if micVp.ispureneumann
%       xi1 = bsxfun(@minus,xi1, get_average(micMesh, micFemspace, xi1));
%     end
%     xi1 = reshape(xi1, dim*ND, 1);
    xi1 = xi1/sqrt(xi1'*W*xi1);
    if isnan(xi1)
        xi1 = zeros(length(xi1),1);
    end
    RBspace = xi1;
end

% Add other basis functions using greedy procedure
for k = 2:NRB
    disp(k)
    NRB_count = size(RBspace,2);
    for q = 1:Q
        for i = 1:NRB_count
            for j = 1:NRB_count
                AqRB(i,j,q) = RBspace(:,i)'*Aq{q}*RBspace(:,j);
            end
            FqRB(i,:,q) = Fq{q}'*RBspace(:,i);
        end
    end
    xik = zeros(k-1, length(SetTrain));
    aposteriori_error = zeros(NTrain, 1);
    res_RB = zeros(NTrain, 1);
    relres_RB = zeros(NTrain, 1);
    relative_error = zeros(NTrain, 1);
    for i = SetTrain
        disp(i)
        ARB = zeros(NRB_count);
        FRB = zeros(NRB_count, micVp.NF);
        for q = 1:Q
            ARB = ARB + Theta(q,i)*AqRB(:,:,q);
            FRB = FRB + Theta(q,i)*FqRB(:,:,q);
        end
        xik(1:k-1, i) = ARB\FRB(:,ChiTrain{2}(i));
        if isnan(xik(:,i))
            xik(:,i) = zeros(k-1,1);
        end
        A = sparse(dim*ND, dim*ND);
        F = sparse(dim*ND, 1);
        for q = 1:Q
            A = A + Theta(q,i)*Aq{q};
            F = F + Theta(q,i)*Fq{q}(:,ChiTrain{2}(i));
        end
        residual = A*RBspace*xik(:,i) - F;
        res_RB(i) = norm(residual);
        relres_RB(i) = norm(residual)/norm(F);
        aposteriori_error(i) = sqrt(residual'*(W\residual));
        relative_error(i) = aposteriori_error(i)/sqrt(F'*(W\F));
        if isnan(aposteriori_error(i))
            aposteriori_error(i) = Inf;
        end
    end
%     [e_new, i_new] = max(aposteriori_error);
    [e_new, i_new] = max(relative_error);
    RBerror(k-1) = e_new;
%     disp(e_new^2);
%     if e_new^2 < tolRB
%         break;
    disp(e_new);
    if e_new < tolRB
        break;
    else
        A = sparse(dim*ND,dim*ND);
        F = sparse(dim*ND,1);
        for q = 1:Q
            A = A + Theta(q,i_new)*Aq{q};
            F = F + Theta(q,i_new)*Fq{q}(:,ChiTrain{2}(i_new));
        end
        xi_new = A\F;
%         xi_new = reshape(xi_new, ND, dim);
%         if micVp.ispureneumann
%             xi_new = bsxfun(@minus,xi_new, get_average(micMesh, micFemspace, xi_new));
%         end
%         xi_new = reshape(xi_new, dim*ND, 1);
        distance = 0;
        for l = 1:k-1
            distance = distance + xi_new'*W*RBspace(:,l)*RBspace(:,l);
        end
        RL = xi_new - distance;
    end
    new_basis_function = RL/sqrt(RL'*W*RL);
    if ~isnan(new_basis_function)
        RBspace = [RBspace, new_basis_function];
    end
    SetTrain = setdiff(SetTrain, i_new);
end    
end

