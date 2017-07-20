%%%%%%%%%%%

%%%%%%%%%%%
% y= [0.0038 0.0080 0.0140 0.0256 0.0522 0.0742 0.1114 0.1418 0.2094 0.3596]';
% y= [0.0118 0.0395 0.1265 0.2532 0.5690];2
% y=a;
% a
clear s;

% x=[1:25,30,35,50];
na=[2 3 5 8 10 12 15 18 20 25 30 50 100];
clear ress;
ress=cell(length(na),1);
pp=[];
p=.3;
% for p=[.35]
   for p=[0.001 .01 .1 .2 .3 .4 .005 .05 .5 .6 .7 .8 .9 .0025 .025 .35 .45 .075 .0075] 
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
        %%%%%%%GA unequal%%%%%%%%%
        clear y;
        y=InterpolateGAProbs(na(i));
        %%%%%%%%%%%%%%%%%%%%%%%
        
        %%%bestofgen%%%%%%%%%%%%%
        %     y=bestofgen{end};
        %     y=y';
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%equal%%%%%%%%%%%%%
%                 y=ones(1,na(i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %        p=0.6;
%         432/2,2,0,1
        ress{i}=CA_Functions2(y,length(y),432,2,1,1,300,p,10);  %probs,numants,numits*10000,width,infEnergy
        %     yy=ress(i).atFace(end,:);
        % pts('size=',i,' (G1,G2,G3)=(',G1,',',G2(l),',',G3(l),')');
        % save([pwd,'\pauseVarypell=100\pause',num2str(i),'.mat'],'res', '-v7.3')
        
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
        
        saveFold = [pwd,'\p',num2str(p),feq,'pell200'];
        if(~exist(saveFold))
            mkdir(saveFold);
        end
        save([saveFold,'\res',num2str(res.numants),'.mat'],'res', '-v7.3')
    end
end
beep;
%% plot gini vs. number of ants NOT UPDATED WITH NEW GINI CALC NUMBER
figure(12); hold on;
titl='Gini vs. Number of ants (48 h)';
xlabel('N (ants)'); ylabel('Gini');
eq=1; %% zero for inequal work
p='0.4';
ptt=[.001 .0025 .005 .0075 .01 .025 .05 .075 .1 .2 .3 .4 .5];
ptt=sort(ptt,'ascend');
outer=jet(length(ptt));
colormap jet
for z=1:length(ptt)
    p=ptt(z);
    pp=num2str(p);
    for eq=[0] %% zero for inequal work
         clear res;
        xText=60;
        if(eq)
            
            fold = ['\p',pp,'eqpell200\'];
            %     fold ='\eqResSavesInfEnergyPell600\';
            %     fold ='\eqResSavesRegEnergyPell200\';
            textz = 'Equal';
            yText = .25;
            mark = '-';
            co = [1 0 0];
        else
            fold =['\p',pp,'uneqpell200\'];
            %     fold ='\uneqResSavesInfEnergyPell600\';
            %     fold ='\uneqResSavesRegEnergyPell200\';
            mark = '-';
            yText = .54;
            textz = 'Unequal';
            co =[184 77 157]./255;
        end
        flist=dir([pwd,fold,'*.mat']);
        for i = 1:length(flist);
            num(i)=sscanf(flist(i).name,'res%d');
        end
        [~,order] = sort(num);
        flist(:)=flist(order);
        
        num=zeros(1,length(flist));
        G=zeros(1,length(flist));
        for i = 1:length(flist);
            load([pwd,fold,flist(i).name]);
            num(i)=sscanf(flist(i).name,'res%d');
            pells= sum(res.markMatr(:,2:end)); %25=pause2dig
            pells2=sum(pells); %cm/timeframe
            %             pells2=pells2/((stopTime-startTime)); % cm/step
            %             pells2=pells2/dt; %(cm/step) /(sec/step)= cm/sec
            G(i)=Gini(pells);
        end
        
        if(res.infEnergy)
            En=' E=\infty';
        else
            En=[];
        end
%         title([titl,En]);
        
        plot(num,G,mark,'color',outer(z,:),'markersize',8,'markeredgecolor',co,'linewidth',2.5);
        p
        
        
    end
%     pause
end
set(gca,'box','on','linewidth',2,'fontsize',25);
%% plot growth rate vs ants with both modes (*PAPER QUALITY*) tunneltip=10?
clear all;

figure(11); hold on;
% titl= 'Tunnel Growth Rate vs. N';
xlabel('{\itn}'); ylabel('{\itV} (cm/h)');
% p='0.3'
ptt=[.001 .0025 .005 .0075 .01 .025 .05 .075];
outer=jet(length(ptt));
for z=1:length(ptt)
    p=ptt(z);
    pp=num2str(p);
    for eq=[0 1] %% zero for inequal work
        clear res;
        xText=60;
        if(eq)
            
            fold = ['\p',pp,'eqpell200\'];
            %     fold ='\eqResSavesInfEnergyPell600\';
            %     fold ='\eqResSavesRegEnergyPell200\';
            textz = 'Equal';
            yText = .25;
            mark = 's-';
            co = [1 0 0];
        else
            fold =['\p',pp,'uneqpell200\'];
            %     fold ='\uneqResSavesInfEnergyPell600\';
            %     fold ='\uneqResSavesRegEnergyPell200\';
            mark = 'o-';
            yText = .54;
            textz = 'Unequal';
            co =[184 77 157]./255;
        end
        flist=dir([pwd,fold,'*.mat']);
        for i = 1:length(flist);
            num(i)=sscanf(flist(i).name,'res%d');
        end
        [~,order] = sort(num);
        flist(:)=flist(order);
        
        antNum=zeros(1,length(flist));
        V=zeros(1,length(flist));
        startTime=.25;
        stopTime = 3+startTime; %in hours
        startTime=startTime*7200;
        stopTime = stopTime*7200;
        
        dt=.5;
        %width*pellet2grow=tunnel tip growth of .5 cm
        %width*pellet2grow*2 = pells/cm
        for i =1:length(flist)
            flist(i).name;
            load([pwd,fold,flist(i).name]);
            antNum(i)=res.numants;
            pellet2grow=res.pell2grow;
            pellPerCm = pellet2grow*2*2;
            %         pells= floor(cumsum(res.markMatr(startTime:stopTime,2:end))/res.pause2dig)/pellPerCm; %25=pause2dig
            %         pellFinal = sum(pells(end,:));
            %         pells2= sum(pells,2);
            %         V(i)= mean(diff(pells2).*3600/(dt));
            %
            pells= sum(res.markMatr(startTime:stopTime,2:end))/res.pause2dig/pellPerCm; %25=pause2dig
            pells2=sum(pells); %cm/timeframe
            pells2=pells2/((stopTime-startTime)); % cm/step
            pells2=pells2/dt; %(cm/step) /(sec/step)= cm/sec
            V(i)= pells2.*3600; %cm/sec*3600s/h
            
            %     pause
        end
        if(res.infEnergy)
            En=' E=\infty';
        else
            En=[];
        end
        % title([titl,' (',num2str((stopTime-startTime)/7200),'h)',En]);
        %%%UNCOMMENT%%%%
        %         plot(antNum,V,mark,'LineWidth',2,'color',co,'markersize',8,'linewidth',4);
        plot(antNum,V,mark,'color',outer(z,:),'markersize',8,'markeredgecolor',co,'linewidth',3);
        
        if (z==1)
            text(xText,yText,textz,'fontsize',35,'color',co); %%%UNCOMMENT%%%%
        end
        % set(xlabh,'Position',get(xlabh,'Position')+[0 .01 0]);
    end
end
% x=2.6 y = 1
set(gcf,'Position',[100,100, 1300/1.2, 500/1.2])

xlabh=get(gca,'xlabel');
axis([0 100 0 0.7]);
set(gca,'box','on','linewidth',2,'xtick',[0 50 100],'ytick',[.3 .6],'fontsize',35);
%% plot q_bar vs rho_bar ****PAPER QUALITY*************SET tunneltip=10
figure(14); hold on;
titl='Traffic (CA)';
% p='0.35'
% h1=xlabel('$\bar{\rho}$ (ants/site)');
% h2=ylabel('$\bar{q}$ (ants/min)');
fsize=24;
h1=xlabel('$\bar{\rho}$ (ants/BL)');
% h2=ylabel('$\bar{q}$ ($\frac{ants}{BL{\cdot}min}$)');
h2=ylabel('$\bar{q}$ (ants/(BL$\cdot$min))','fontsize',fsize);

set(h1,'interpreter','Latex','FontWeight','bold');
set(h2,'interpreter','Latex','FontWeight','bold','fontsize',fsize);
% eq=0; %% zero for inequal work
% ptt=[0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075];
ptt=[0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.2 0.3 ];
ptt=sort(ptt,'descend');
for z=1:length(ptt)
    p=ptt(z);
    pp=num2str(p);
    for eq = [0]
        clear res;
        if(eq)
            %     fold ='\eqResSavesInfEnergy\';
            %     fold ='\eqResSavesInfEnergyPell600\';
            fold = ['\p',pp,'eqpell200\'];
            mark = 's';
            lw=2;
        else
            %     fold ='\uneqResSavesInfEnergy\';
            %     fold ='\uneqResSavesInfEnergyPell600\';
            fold =['\p',pp,'uneqpell200\'];
            mark = 'o';
            lw = 1.9;
        end
        
        
        flist=dir([pwd,fold,'*.mat']);
        for i = 1:length(flist);
            num(i)=sscanf(flist(i).name,'res%d');
        end
        [~,order] = sort(num);
        flist(:)=flist(order);
        ants=[2 3 5 8 10 12 15 18 20 25 30 50 100];
        numcolors = length(ants*3);
        
        stopTime =2.75;
        startTime=.25;
        minz = 60*(stopTime-startTime); %%%uncomment
        %   minz = 60*(stopTime-startTime);
        stopTime = stopTime*7200;
        startTime=startTime*7200;
        
        outer= jet(length(ptt));
        cmap = fire(numcolors);
        for i =1:length(flist);
            flist(i).name;
            
            load([pwd,fold,flist(i).name]);
            %         singleLane = sum(res.occupied(startTime+1:stopTime,:),2)./2; %/2 for 2 lanes
            singleLane = sum(res.density(startTime+1:stopTime,:),2)./2; %/2 for 2 lanes
            singleFlow = sum(res.flow(startTime+1:stopTime,:),2)./2;
            tunLen=(res.tunLength(startTime+1:stopTime,2));
            ts=120*minz; %timesteps to average over 120ts/min *60min/hour
            %     ss=size(res.occupied(1:stopTime,1));
            
            %         newq = reshape(singleFlow(1:end),ts,size(singleFlow(1:end),1)/ts);
            %         tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
            %         newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts);
            
            
            a=sum(sum(res.markMatr(startTime:stopTime,2:end))/res.pause2dig);
            newq=a/ts;
            tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
            newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts)./tunLen;
            
            
            
            qm = mean(newq);
            tunLen=mean(tunLen);
            rhom = mean(newrho);
            
            %             rhom = mean(rhom);
            %             qm = mean(qm);
            %           rhom=mean(singleLane./tunLen);
            %           qm=mean(singleFlow);
            %%%%%%%%%%%%%%%%%%%%WHAT SHOULD WORK?%%%%%%%%%%%%%%%%%%%%%%%
            %             rhom = mean(singleLane).*tunLen(1)./tunLen;
            %             qm = mean(singleFlow);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plot(rhom,qm,mark,'MarkerSize',15,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',outer(z,:),'LineWidth',lw);
            %             plot(rhom,qm,mark,'MarkerSize',15,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k','LineWidth',lw);
            %                     plot(rhom,qm,mark,'MarkerSize',6,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',outer(str2double(p)*10,:),'LineWidth',1);
            
        end
        haxis = gca;
        set(gca,'fontsize',55);
        figText(gcf,16)
        colormap(cmap);
        
        caxis([1 numcolors])
        
        cbarHandle = colorbar('YTick',...
            1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors,...
            'YTickLabel',ants, 'YLim', [1 numcolors]); %[]=ants
        
        % set(gca,'XTick',[0.25 0.5],'YTick',[0 .1 .2],'XLim',[0.045 0.55],'YLim',[0,0.2]);
        set(gca,'XTick',[0 0.15 0.3],'YTick',[ .1 .2],'fontsize',fsize);
        set(cbarHandle,'fontsize',18');
        axis([0 .5 0 .1])
        if(res.infEnergy)
            En=' E=\infty';
        else
            En=[];
        end
        % title([titl,' (',num2str(stopTime/7200),'h)',En]);
        % title([titl,En]);
        set(gca,'box','on','linewidth',2);
    end
    p
    pause;
    
    
end
%%
%% plot q_bar vs rho and reversal
figure(15); hold on;
% titl='Traffic (CA)';
% p='0.35'
% h1=xlabel('$\bar{\rho}$ (ants/site)');
% h2=ylabel('$\bar{q}$ (ants/min)');
fsize=24;
h1=xlabel('$\bar{\rho}$ (ants/BL)');
% h2=ylabel('$\bar{q}$ ($\frac{ants}{BL{\cdot}min}$)');
h2=ylabel('$\bar{q}$ (ants/(BL$\cdot$min))','fontsize',fsize);

set(h1,'interpreter','Latex','FontWeight','bold');
set(h2,'interpreter','Latex','FontWeight','bold','fontsize',fsize);
% eq=0; %% zero for inequal work
% ptt=[0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075];
ptt=[0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.2 0.35 ];
ptt=sort(ptt,'descend');
for z=1:length(ptt)
    p=ptt(z);
    pp=num2str(p);
    for eq = [0]
        clear res;
        if(eq)
            %     fold ='\eqResSavesInfEnergy\';
            %     fold ='\eqResSavesInfEnergyPell600\';
            fold = ['\p',pp,'eqpell200\'];
            mark = 's';
            lw=2;
        else
            %     fold ='\uneqResSavesInfEnergy\';
            %     fold ='\uneqResSavesInfEnergyPell600\';
            fold =['\p',pp,'uneqpell200\'];
            mark = 'o';
            lw = 1.9;
        end
        
        
        flist=dir([pwd,fold,'*.mat']);
        for i = 1:length(flist);
            num(i)=sscanf(flist(i).name,'res%d');
        end
        [~,order] = sort(num);
        flist(:)=flist(order);
        ants=[2 3 5 8 10 12 15 18 20 25 30 50 100];
        numcolors = length(ants*3);
        
        stopTime =2.75;
        startTime=.25;
        minz = 60*(stopTime-startTime); %%%uncomment
        %   minz = 60*(stopTime-startTime);
        stopTime = stopTime*7200;
        startTime=startTime*7200;
        
        outer= jet(length(ptt));
        cmap = fire(numcolors);
        for i =1:length(flist);
            flist(i).name;
            
            load([pwd,fold,flist(i).name]);
            %         singleLane = sum(res.occupied(startTime+1:stopTime,:),2)./2; %/2 for 2 lanes
            singleLane = sum(res.density(startTime+1:stopTime,:),2)./2; %/2 for 2 lanes
            singleFlow = sum(res.flow(startTime+1:stopTime,:),2)./2;
            tunLen=(res.tunLength(startTime+1:stopTime,2));
            ts=120*minz; %timesteps to average over 120ts/min *60min/hour
            %     ss=size(res.occupied(1:stopTime,1));
            
            %         newq = reshape(singleFlow(1:end),ts,size(singleFlow(1:end),1)/ts);
            %         tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
            %         newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts);
            
            
            a=sum(sum(res.markMatr(startTime:stopTime,2:end))/res.pause2dig);
            newq=a/ts;
            tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
            newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts)./tunLen;
            
            
            
            qm = mean(newq);
            tunLen=mean(tunLen);
            rhom = mean(newrho);
            
            %             rhom = mean(rhom);
            %             qm = mean(qm);
            %           rhom=mean(singleLane./tunLen);
            %           qm=mean(singleFlow);
            %%%%%%%%%%%%%%%%%%%%WHAT SHOULD WORK?%%%%%%%%%%%%%%%%%%%%%%%
            %             rhom = mean(singleLane).*tunLen(1)./tunLen;
            %             qm = mean(singleFlow);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plot(rhom,qm,mark,'MarkerSize',15,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',outer(z,:),'LineWidth',lw);
            %             plot(rhom,qm,mark,'MarkerSize',15,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k','LineWidth',lw);
            %                     plot(rhom,qm,mark,'MarkerSize',6,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',outer(str2double(p)*10,:),'LineWidth',1);
            
        end
        haxis = gca;
        set(gca,'fontsize',55);
        figText(gcf,16)
        colormap(cmap);
        
        caxis([1 numcolors])
        
        cbarHandle = colorbar('YTick',...
            1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors,...
            'YTickLabel',ants, 'YLim', [1 numcolors]); %[]=ants
        
        % set(gca,'XTick',[0.25 0.5],'YTick',[0 .1 .2],'XLim',[0.045 0.55],'YLim',[0,0.2]);
        set(gca,'XTick',[0 0.15 0.3],'YTick',[ .1 .2],'fontsize',fsize);
        set(cbarHandle,'fontsize',18');
        axis([0 .5 0 .1])
        if(res.infEnergy)
            En=' E=\infty';
        else
            En=[];
        end
        % title([titl,' (',num2str(stopTime/7200),'h)',En]);
        % title([titl,En]);
        set(gca,'box','on','linewidth',2);
    end
    p
    pause;
    
    
end