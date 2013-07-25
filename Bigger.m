function result = Bigger(smallMatrix,increaseFactor)

[m,n] = size(smallMatrix);

temp = [];
result = [];
for i = 1:m
    for j = 1:increaseFactor
        temp = [temp; smallMatrix(i,:)];
    end
end

result = temp;
% for i = 1:n
%     for j = 1:increaseFactor
%         result = [result temp(:,i)];
%     end
% end

