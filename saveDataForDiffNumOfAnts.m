clear s;
[filez,fold]=uigetfile(fullfile('D:\Projects\Ant_CA_GA\results','*.mat'));
na=[2 3 5 8 10 12 15 18 20 25 30 50 100];
clear ress;
ress=cell(length(na),1);
% filez=uigetfile(fullfile('D:\Projects\Ant_CA_GA\results','*.mat'));
load(fullfile(fold,filez));
for p=[prob2turn]
%    for p=[0.001 .01 .1 .2 .3 .4 .005 .05 .5 .6 .7 .8 .9 .0025 .025 .35 .45 .075 .0075] 
    clear res ress;
    ress=cell(length(na),1);
    for i=1:length(na)
        %%%%%%%%unequal%%%%%%%%%
%         clear s
%         s.x =1;
%         s= SetProbabilities(pwd,na(i),s,1,0); %9 - probabilities!!! from averaged Lorenz curve for 12 hours
%         y=s.prob;%/sum(s.prob.^(1.75));
%         % %         y=y./sum(y);
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%GA unequal%%%%%%%%%
%         clear y;
%         pp=InterpolateGAProbsFromProb(na(i),bestofgen{end});
        %%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%equal%%%%%%%%%%%%%
        pp=ones(1,na(i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ress{i}=CA_FunctionsWill(pp,length(pp),432,TW,energyMult,1,rechargeSteps,prob2turn,tuntip);
%         ress{i}=CA_Functions2(y,length(y),432,2,1,1,300,p,10);  %probs,numants,numits*10000,width,infEnergy
% %         %     yy=ress(i).atFace(end,:);
% %         % pts('size=',i,' (G1,G2,G3)=(',G1,',',G2(l),',',G3(l),')');
% %         % save([pwd,'\pauseVarypell=100\pause',num2str(i),'.mat'],'res', '-v7.3')
     i   
    end
    feq=0;
    a=all(ress{1}.prob(1,1)==ress{1}.prob(1,:));
    if all(ress{1}.prob(1,1)==ress{1}.prob(1,:))
        feq='eq';
    else
        feq='uneq';
    end
    feq
    for i=1:length(na)
        res=ress{i};
        
        saveFold = ['D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\diffN',feq];
        if(~exist(saveFold))
            mkdir(saveFold);
        end
        save([saveFold,'\res',num2str(res.numants),'.mat'],'res', 'pp', '-v7.3')
    end
end
beep;
