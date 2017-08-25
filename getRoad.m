%%%%%%%regular%%%%%%%%
nvars = 30;
%25 30 35 40 45 50 55 60 65 70
%%%%%%%%%%%%%%%%%%%%%%
totgens = 10;
popsize = 7*nvars;
tuntip=10;
prob2turn = .3;
rechargeSteps = 10;
numIts=432;
energyMult=1;
TW=2;
% prob=sort(bestofgenOUT{end}/sum(bestofgenOUT{end}));
% [res,road]=CA_FunctionsWill(prob,length(prob),numIts,TW,...
%     energyMult,1,rechargeSteps,prob2turn,tuntip);
% filez=uigetfile('D:\Projects\Ant_CA_GA\results');
filez='D:\Projects\Ant_CA_GA\results\unequal seed 5-100 tw=2 gen=20\N=30_tw=2_2017-08-17-00-03.mat';
load(filez);
[res,road]=CA_FunctionsWill(bestofgen{end},length(bestofgen{end}),numIts,TW,...
    energyMult,1,rechargeSteps,prob2turn,tuntip);
% 
% cmap=brewermap(32*2,'RdGy');
% colormap(cmap);
figure(12);
cmap=brewermap(32*2,'RdGy');
colormap(cmap);
sIx=1999; %start idx for plotting
subplot(3,1,1);
imagesc(road(:,:,sIx));
subplot(3,1,2);
imagesc(road(:,:,sIx+1));
subplot(3,1,3);
imagesc(road(:,:,sIx+2));