for(jj=[0.01:0.01:0.1,.1:.025:1])
%%%%%%%regular%%%%%%%%
nvars = 30;
%25 30 35 40 45 50 55 60 65 70
%%%%%%%%%%%%%%%%%%%%%%
totgens = 10;
popsize = 7*nvars;
tuntip=10;
% prob2turn = .06;
rechargeSteps = 10;
numIts=432;
energyMult=1;
TW=2;
typez=1; %0=equal,1=unequal
% [res,road]=CA_FunctionsWill(prob,length(prob),numIts,TW,...
%     energyMult,1,rechargeSteps,prob2turn,tuntip);
% filez=uigetfile('D:\Projects\Ant_CA_GA\results');
filez='D:\Projects\Ant_CA_GA\results\unequal seed 5-100 tw=2 gen=20\N=30_tw=2_2017-08-17-00-03.mat';
load(filez);
prob2turn = jj;
if(typez==0)
pp=ones(nvars,1);
else
pp=sort(bestofgenOUT{end}/sum(bestofgenOUT{end}));
end
[res,roadFull]=CA_FunctionsWill(pp,length(pp),numIts,TW,...
    energyMult,1,rechargeSteps,prob2turn,tuntip);
%% save out data
dateFormat='dd-mm-yyyy_HHMM_';
a=datestr(now,dateFormat);
roadSmall=roadFull(:,:,1000:1100);
% fnameOut=fullfile(pwd,[a,'roadData.mat']);
fnameOut=fullfile(pwd,'results','bahni2',['type_',num2str(typez),'_R_',num2str(prob2turn),'.mat']);
save(fnameOut,'roadFull','roadSmall','res');

end

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

