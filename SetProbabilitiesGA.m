function prob = SetProbabilitiesGA(numants)
%SETPROBABILITIES Summary of this function goes here
%   Detailed explanation goes here
%     cd(workfolder);
%     load('InterpLorenzC.mat');
%     load('mode3.mat','pop')

%     x1=[0 (1:4)/4];
%     y1=[0 pop((end-3):end)];%percent
%     y1=y1/sum(y1);
%     y1=cumsum(y1);

    load('InterpLorenzC.mat');
    x1(:,1)=xq;%prob
    y1(:,1)=avgXstd(:,1)';%percent
    
    oneant=1/numants;
    for i=1:numants
        %p=interp1(x1,y1,i*oneant)-interp1(x1,y1,(i-1)*oneant);
        p=interp1(x1,y1,i*oneant)-interp1(x1,y1,(i-1)*oneant);
%         systemState.prob(i,1)=p*coeffP;
        prob(i,1)=p;
    end

end


