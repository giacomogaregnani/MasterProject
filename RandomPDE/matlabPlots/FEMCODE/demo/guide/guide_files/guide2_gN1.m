function neumann = guide2_gN1(x, normals)
neumann = cos(2*pi*x(:,1)) .* cos(2*pi*x(:,2)) .* normals(:,1) - ...
    sin(2*pi*x(:,1)) .* sin(2*pi*x(:,2)) .* normals(:,2);
end