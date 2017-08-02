function [ state ] = gaplotlorenzcurve( options,state,flag,bestofgen,numIts,rec,ptt,tuntip,tw,emult )%NORMAL
%gaplotlorenzcurve plotting function for GA connected to darias CA model
% function [ state ] = gaplotlorenzcurve( options,state,flag,bestofgen,numIts)%%WITH rech and ptt




a=evalin('base','bestofgen');
ni=evalin('base','numIts');
na=evalin('base','nvars');
% prob = a{state.Generation+1};

prob = a{state.Generation+1}(1:na); %30 for 30 ants

%%%%%%%%%%rec and ptt%%%%%%%%%%%
% rec = a{state.Generation+1}(31);
% ptt = a{state.Generation+1}(32);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=CA_FunctionsWill(prob,length(prob),ni,tw,emult,1,rec,ptt,tuntip);
f=sum(r.markMatr(:,2:end));
[ginM,gxy]=Gini(f);
plot(gxy(:,1),gxy(:,2),'.-','linewidth',2,'markersize',10);
hold on;
load('antfigData10.mat');
set(gcf,'renderer','openGL');

cc=parula(13);
ccl=cc(2,:);
plot(antfigx,antfigy,'-','color',ccl,'linewidth',2,'markersize',20);
patch(antfigxP,antfigyP,ccl,'facealpha',.15,'edgecolor','none');
load('antfigData1.mat');
plot(antfigx,antfigy,'r-','linewidth',2);
patch(antfigxP,antfigyP,'r','facealpha',.15,'edgecolor','none')
axis([0 1 0 1]);
text(.2,.8,{['{\itGen} = ',num2str(state.Generation)];...
    ['G=',num2str(ginM,3)];['R=',num2str(rec)];['P=',num2str(ptt)]},'Color','k','FontSize',18)
set(gca,'xscale','log','yscale','log');
axis tight;
hold off;
end

