load('mGoodStableRule1PingRej-split_by_Day_BehResp_and_Area.mat')
%create an ERP for a specific day using the session_erp func
session_erp(monkey,1,7);

%add day to the title of the figure in the erp