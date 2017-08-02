bestofgen = bestofgen(~cellfun(@isempty,bestofgen));
bestofgenOUT=cell(1,length(bestofgen));
% y=bestofgen{end};
clear res;

nvars = 100;
%%%%%%%%%%%%%%%%%%%%%%
totgens = 20;
popsize = 3;
tuntip=10;
prob2turn = .3;
rechargeSteps = 10;
numIts=432;
energyMult=1;
TW=5;
parfor i=1:length(bestofgen)
    y=bestofgen{i};
    
    % res(i) =CA_FunctionsWill(y,length(y),numIts,2,0);  %probs,numants,numits*10000,width,infEnergy
    % f=sum(res(i).markMatr(:,2:end));
    res =CA_FunctionsWill(y,length(y),numIts,TW,energyMult,1,rechargeSteps,prob2turn,tuntip);  %probs,numants,numits*10000,width,infEnergy
    
    
    f=sum(res.markMatr(:,2:end));
    bestofgenOUT{i}=f;
end

%% plot how bestofgen gini (expected) changes
if(exist('bestofgen',1))
    figure(8); hold on;
    title('GA Gini (expected) vs. Generations (6h)');
    xlabel('Generations');
    ylabel('Gini');
    gin=zeros(1,length(bestofgen));
    gen=1:length(bestofgen);
    for i = 1:length(bestofgen)
        [gin(i),~]=Gini(bestofgen{i});
    end
    plot(gen,gin,'.-');
end

%% plot how bestofgen gini (measured) changes
if(exist('bestofgenOUT',1))
    figure(8); hold on;
    title('GA Gini (measured) vs. Generations (6h)');
    xlabel('Generations');
    ylabel('Gini');
    gin=zeros(1,length(bestofgenOUT));
    gen=1:length(bestofgenOUT);
    for i = 1:length(bestofgenOUT)
        [gin(i),~]=Gini(bestofgenOUT{i});
    end
    plot(gen,gin,'.-');
end

%% plot how bestofgen gini changes
if(exist('bestofgenOUT',1) && exist('bestofgen',1) )
    figure(8); hold on;
    title('GA Gini vs. Generations (6h)');
    xlabel('Generations');
    ylabel('Gini');
    gin=zeros(1,length(bestofgenOUT));
    gen=1:length(bestofgenOUT);
    clear ginE ginM
    for i = 1:length(bestofgenOUT)
        [ginE(i),~]=Gini(bestofgen{i});
        [ginM(i),~]=Gini(bestofgenOUT{i});
    end
    plot(gen,ginE,'.-');
    plot(gen,ginM,'.-');
    set(gca,'xlim',[1,length(ginM)])
    legend({'expected','measured'});
end

%% plot how bestofgen gini changes movie
ftext=[];

mov=0;
last=0;
if mov
    fold = [uigetdir,'\'];
end

if(exist('bestofgenOUT',1) && exist('bestofgen',1) )
    figure(1);
    title('Lorenz Curve (6h)');
    xlabel('Generations');
    ylabel('Gini');
    clear gxy ginM
    gin=zeros(1,length(bestofgenOUT));
    gen=1:length(bestofgenOUT);
    if ~last
        hold on;
    end
    for i = 1:length(bestofgenOUT)
        [ginM(i),gxy{i}]=Gini(bestofgenOUT{i});
    end
    % plot(gen,ginE,'.-');
    % plot(gen,ginM,'.-');
    % legend({'expected','measured'});
end

for i = 1:length(gxy)
    plot(gxy{i}(:,1),gxy{i}(:,2),'o-');
    title('Lorenz Curve (6h)');
    
    xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}')
    text(.2,.8,['{\itGen} = ',num2str(i)],'Color','k','FontSize',18)
    
    if mov
        saveas(gcf,[fold, num2str(i)],'png')
    end
    % ffmpeg -framerate 4 -i %d.png -c:v libx264 -r 30 out.mp4
    pause(.1)
    % !cd A:\GA_equal\
    % ffmpeg -framerate 4 -i %d.png -c:v libx264 -r 30 out.mp4
end
if mov
    system(['ffmpeg -framerate 4 -i ',fold,'%d.png -c:v mpeg4 -vtag xvid -r 30 -q 0',fold,'out.avi']);
end

%% plot lorenz over ant curve

[GG,gxy]=Gini(bestofgenOUT{end});
figure(2); hold on;
xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}')
plot(gxy(:,1),gxy(:,2),'-o','linewidth',2);
set(gca,'XScale','linear','yscale','linear','xlim',[0.01 1],'ylim',[0.0001 1]);
GG
%  ax = gca;

%   ax.XScale='log';
%   ax.YScale='log';
% mut=.25 p=100 g=50 G=.635
%%
set(gca,'XScale','log','yscale','log','xlim',[0.01 1],'ylim',[0.0001 1]);
axis tight
%%
set(gca,'XScale','linear','yscale','linear','xlim',[0.01 1],'ylim',[0.0001 1]);
%% PLOT ANT DATA

%load ant data
% h=openfig('antLorenz')

set(gcf,'renderer','openGL');
cc=parula(13);
ccl=cc(2,:);
hold on;
load('A:\GA\CA model v2\antfigData10.mat');
plot(antfigx,antfigy,'-','color',ccl,'linewidth',2);
patch(antfigxP,antfigyP,ccl,'facealpha',.15,'edgecolor','none');
load('A:\GA\CA model v2\antfigData1.mat');
plot(antfigx,antfigy,'r-','linewidth',2);
patch(antfigxP,antfigyP,'r','facealpha',.15,'edgecolor','none')
h=gcf;
legend('show');
drawnow
set(h.Children(2).Children(1).Annotation.LegendInformation,'IconDisplayStyle','off')
set(h.Children(2).Children(3).Annotation.LegendInformation,'IconDisplayStyle','off')
% legend('');
xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}')
legend('Ant Data .1', 'Ant Data .01');
%% load
load('A:\energydiv4rechargesteps300\GAdat');


%  load('A:\randg33p100m15P2D1\GAdat');
%  load('A:\oldGaRes\GA_unequal6h\GAdat');
%  load('A:\oldGaRes\GA_equal6h\GAdat');

%% bar plot of
fold= [uigetdir,'\'];
for i = 1:length(bestofgenOUT)
    
    bar(1:30,sort(bestofgenOUT{i})/sum(bestofgenOUT{i}));
    axis([0,31 0 .3]);
    xlabel('ant');
    ylabel('work probability');
    title(['Work Probabilities Gen=',num2str(i)]);
    
    saveas(gcf,[fold, num2str(i)],'png');
    
    pause(.05);
end
system(['ffmpeg -framerate 4 -i "',fold,'%d.png" -c:v libx264 -r 30 ',fold,'out.mp4']);
%%
