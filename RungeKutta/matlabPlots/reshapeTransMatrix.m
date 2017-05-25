function TransMatrix = reshapeTransMatrix(TransMatrix) 

if size(TransMatrix, 1) ~= size(TransMatrix, 2)
    [~, minDim] = min(size(TransMatrix));
    if minDim == 1
       nRowsFill = size(TransMatrix, 2) - size(TransMatrix, 1);
       TransMatrix = [TransMatrix; zeros(nRowsFill, size(TransMatrix, 2))];
    else
       nColFill = size(TransMatrix, 1) - size(TransMatrix, 2);
       TransMatrix = [TransMatrix, zeros(nColFill, size(TransMatrix, 1))'];
    end
end