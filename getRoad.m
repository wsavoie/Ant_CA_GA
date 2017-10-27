% jj=[0.01:0.01:0.1,.1:.025:1];
count = 0;
for(tt=[0 1])
    for(jj=[0.01:0.01:0.1,.1:.025:1])
        clearvars -except jj count tt
        count=count+1;
        %%%%%%%regular%%%%%%%%
        nvars = 30;
        %25 30 35 40 45 50 55 60 65 70
        %%%%%%%%%%%%%%%%%%%%%%
        
        typez=tt; %0=equal,1=unequal
        % [res,road]=CA_FunctionsWill(prob,length(prob),numIts,TW,...
        %     energyMult,1,rechargeSteps,prob2turn,tuntip);
        % filez=uigetfile('D:\Projects\Ant_CA_GA\results');
        % filez='D:\Projects\Ant_CA_GA\results\unequal seed 5-100 tw=2 gen=20\N=30_tw=2_2017-08-17-00-03.mat';
        filez='D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\finEng_12h\rand\N=30_tw=2_2017-10-22-23-01.mat';
        load(filez);
        % prob2turn = .06;
        numIts=432*4;
        
        
        prob2turn = jj;
        if(typez==0)
            pp=ones(nvars,1);
        else
            pp=sort(bestofgen{end}/sum(bestofgen{end}));
        end
        [res,roadFull]=CA_FunctionsWill(pp,length(pp),numIts,TW,...
            energyMult,1,rechargeSteps,prob2turn,tuntip);
        
        % end
        %% save out data
        dateFormat='dd-mm-yyyy_HHMM_';
        a=datestr(now,dateFormat);
        roadSmall=roadFull(:,:,1000:1100);
        % fnameOut=fullfile(pwd,[a,'roadData.mat']);
        fnameOut=fullfile('D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\finEng_12h\reversal data 432x4 its',['type_',num2str(typez),'_R_',num2str(prob2turn),'.mat']);
        save(fnameOut,'roadFull','roadSmall','res');
        pts(count,'/',47*2);
    end
end
%
% while i ~= 1
%% plot out
% figure(12);
% cmap=brewermap(32*2,'RdGy');
% colormap(cmap);
% sIx=1999; %start idx for plotting
% subplot(3,1,1);
% imagesc(road(:,:,sIx));
% subplot(3,1,2);
% imagesc(road(:,:,sIx+1));
% subplot(3,1,3);
% imagesc(road(:,:,sIx+2));

