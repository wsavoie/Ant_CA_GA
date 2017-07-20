function systemState = SetProbabilities(workfolder,numants,systemState,coeffP,all1)
%SETPROBABILITIES Summary of this function goes here
%   Detailed explanation goes here
if all1==0
    cd(workfolder);
    load InterpLorenzC.mat avgXstd xq;
%     load InterpLorenzC.mat avgXstd xq
    x1(:,1)=xq;%prob
    y1(:,1)=avgXstd(:,1)';%percent, variable found in InterpLorenzC
%     y1(:,1)=zeros(21,1);
%     y1(end)=1;
    oneant=1/numants;
    for i=1:numants
        %p=interp1(x1,y1,i*oneant)-interp1(x1,y1,(i-1)*oneant);
        p=interp1(x1,y1,i*oneant)-interp1(x1,y1,(i-1)*oneant);
%         systemState.prob(i,1)=p*coeffP;
        systemState.prob(i,1)=p;
    end

else
   systemState.prob(1:numants,1)=1;
end
end


