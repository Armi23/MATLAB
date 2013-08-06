function result = Bigger(smallMatrix,increaseFactor)

[m,~] = size(smallMatrix);

result = [];
for i = 1:m
    for j = 1:increaseFactor
        result = [result; smallMatrix(i,:)]; %#ok<AGROW>
    end
end
