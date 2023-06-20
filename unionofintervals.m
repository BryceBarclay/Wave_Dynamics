function unionset = unionofintervals(intervals)
% intervals is already sorted with the smallest left end point first

i = 1; s = 1;
while(i <= length(intervals(:,1)))
    left = intervals(i,1); right = intervals(i,2);
    j = i+1;
    if(i < length(intervals(:,1)))
    while(intervals(j,1) <= right)
        if(intervals(j,2) > right)
            right = intervals(j,2);
        end
        j = j+1;
        if(j == length(intervals(:,1))+1)
            break;
        end
    end
    end
    unionset(s,:) = [left,right];
    s = s+1;
    i = j;
end


end