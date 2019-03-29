function pde = define_pde()

pde = struct('ThetaField', @ThetaField, 'tensor', @tensor, 'micTensor', @micTensor, 'micTensorQ', @micTensorQ, 'f', @f, 'exactu', @exactu, 'g_D', @g_D, 'g_N', @g_N);
   
    function z = ThetaField(p,q,parameters)
        NP = size(p,1);
        E = parameters(1) + 0.*p(:,1);
        nu = parameters(2) + 0.*p(:,1);
        switch q
            case 1
                z = 1./(1-nu.^2).*E;
            case 2
                z = 1./(1-nu.^2).*nu.*E;
            case 3
                z = 1./(1-nu.^2).*nu.*E;
            case 4
                z = 1./(1-nu.^2).*E;
            case 5
                z = 1./(1+nu).*E.*0.5;
        end
    end


    function z = tensor(p,i,j,parameters)
        if i == 1 && j == 1
            z = ThetaField(p,1,parameters);
        elseif i == 1 && j == 2
            z = ThetaField(p,2,parameters);
        elseif i == 2 && j == 1
            z = ThetaField(p,3,parameters);
        elseif i == 2 && j == 2
            z = ThetaField(p,4,parameters);
        elseif i == 3 && j == 3
            z = ThetaField(p,5,parameters);
        else
            z = zeros(size(p,1),1);
        end
    end
    
%     function z = tensor(p,i,j)
%         if i == 1 && j == 1
%             z = ThetaField(p,1).*500./(5 + 3.5*sin(2*pi*p(:,1)*64));
%         elseif i == 1 && j == 2
%             z = ThetaField(p,2).*35*cos(4*pi*p(:,2)*64).^2+10;
%         elseif i == 2 && j == 1
%             z = ThetaField(p,3).*35*cos(4*pi*p(:,2)*64).^2+10;
%         elseif i == 2 && j == 2
%             z = ThetaField(p,4).*500./(5 + 3.5*cos(2*pi*p(:,1)*64));
%         elseif i == 3 && j == 3
%             z = ThetaField(p,5).*50*(sin(2*pi*p(:,1)*64).*sin(2*pi*p(:,2)*64)+1) + 10;
%         else
%             z = zeros(size(xmic,1),1);
%         end
%     end
    
    function z = micTensor(p,xmic,i,j,parameters)
        NP = size(xmic,1);
        if i == 1 && j == 1
            z = ThetaField(p,1,parameters)*ones(NP,1);
        elseif i == 1 && j == 2
            z = ThetaField(p,2,parameters)*ones(NP,1);
        elseif i == 2 && j == 1
            z = ThetaField(p,3,parameters)*ones(NP,1);
        elseif i == 2 && j == 2
            z = ThetaField(p,4,parameters)*ones(NP,1);
        elseif i == 3 && j == 3
            z = ThetaField(p,5,parameters)*ones(NP,1);
        else
            z = zeros(NP,1);
        end
    end

%     function z = micTensor(p,xmic,i,j)
%         if i == 1 && j == 1
%             z = ThetaField(p,1).*500./(5 + 3.5*sin(2*pi*xmic(:,1)));
%         elseif i == 1 && j == 2
%             z = ThetaField(p,2).*35*cos(4*pi*xmic(:,2)).^2+10;
%         elseif i == 2 && j == 1
%             z = ThetaField(p,3).*35*cos(4*pi*xmic(:,2)).^2+10;
%         elseif i == 2 && j == 2
%             z = ThetaField(p,4).*500./(5 + 3.5*cos(2*pi*xmic(:,1)));
%         elseif i == 3 && j == 3
%             z = ThetaField(p,5).*50*(sin(2*pi*xmic(:,1)).*sin(2*pi*xmic(:,2))+1) + 10;
%         else
%             z = zeros(size(xmic,1),1);
%         end
%     end

    function z = micTensorQ(xmic, i, j, q)
        NP = size(xmic,1);
        switch q
            case 1
                if i == 1 && j == 1
                    z = ones(NP,1);
                else
                    z = zeros(NP,1);
                end
            case 2
                if i == 1 && j == 2
                    z = ones(NP,1);
                else
                    z = zeros(NP,1);
                end
            case 3
                if i == 2 && j == 1
                    z = ones(NP,1);
                    
                else
                    z = zeros(NP,1);
                end
            case 4
                if i == 2 && j == 2
                    z = ones(NP,1);
                else
                    z = zeros(NP,1);
                end
            case 5
                if i == 3 && j == 3
                    z = ones(NP,1);
                else
                    z = zeros(NP,1);
                end
            otherwise
                if i == 1 && j == 1
                    z = ones(NP,1);
                elseif i == 1 && j == 2
                    z = ones(NP,1);
                elseif i == 2 && j == 1
                    z = ones(NP,1);
                elseif i == 2 && j == 2
                    z = ones(NP,1);
                elseif i == 3 && j == 3
                    z = ones(NP,1);
                else
                    z = zeros(NP,1);
                end
        end
    end
            
    function z = f(p)
        z = [0.*p(:,1) + 0.*p(:,2), 0.*p(:,1) + 0.*p(:,2)];
    end

    function z = exactu(p)
        z = [0.*p(:,1) + 0.*p(:,2), 0.*p(:,1) + 0.*p(:,2)];
    end

    function z = g_D(p)
        z = exactu(p);
    end

    function z = g_N(p, normals)
        zx = zeros(size(p,1),2);
        zx = sum(zx.*normals,2);
        zy = zeros(size(p,1),2);
        zy(:,2) = -1 + 0.*p(:,1) + 0.*p(:,2);
        zy = sum(zy.*normals,2);
        z = [zx, zy];
    end

end