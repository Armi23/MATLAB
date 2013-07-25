function matrixresult = SmoothMatrix(data)

[m,n] =  size(data);
matrixresult = zeros(m,n);
for i = 1:n
    matrixresult(:,i) =  SmoothList(data(:,i));
end
