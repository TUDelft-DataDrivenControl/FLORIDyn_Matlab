function YawDataNew = condenseSOWFAYaw(YawData)
%CONDENSESOWFAYAW Reduce yaw data to what is needed
diff = sum(...
    abs(YawData(2:end-1,2:end) - YawData(1:end-2,2:end)) + ...
    abs(YawData(2:end-1,2:end) - YawData(3:end,2:end)),2);
ind_important = [1;find(diff>0)+1;size(YawData,1)];
YawDataNew = YawData(ind_important,:);
end

