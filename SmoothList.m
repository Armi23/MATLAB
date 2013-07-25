function smoothlist = SmoothList(list)

[m,~] = size(list);

% Zeros must be removed because they are outliers. Find their indices so
% that we can break apart the list at those points. 
for i = 1:m
    if list(i)
        index = i;
        break
    end
end

endindex = 0;
for i = index:m
    if ~list(i)
        endindex = i;
        break
    end
end

% If the end of the list does not have zeros, we only have 2 parts to work
% with, otherwise we have 3 (2 lists of zero and 1 of numbers);
if ~(endindex)
    startOfList = list(1:index-1);
    midOfList = list(index:end);
    smoothlist = [startOfList; smooth(midOfList)];
else
    startOfList = list(1:index-1);
    midOfList = list(index:endindex - 1);
    endOfList = list(endindex:end);
    smoothlist = [startOfList; smooth(midOfList); endOfList];    
end