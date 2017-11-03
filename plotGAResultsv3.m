clear all
% close all
% fold=uigetdir('A:\SmarticleRun\')
% load(fullfile(fold,'amoebaData.mat'));

% load('A:\SmarticleRun\Amoeba_newsquare_1_dead\amoebaData.mat');



%************************************************************
%* Fig numbers:
%* 1. gini vs gens for different numbers of ants
%* 2. best of gen gini(last) vs. number of ants
%* 3. number of ants, q vs rho *use 7*
%* 4. (gen(10)-gen(1)) vs. number of ants
%* 5/6. Tunnel length excavated vs. time  *AND* per ant FOR GA vals
%* 7. *fixed* workload distribution q vs rho *use 16* for multires files
%* 8. get cluster dissolution
%* 9. excavation amount for results for different single runs  V vs reversal
%*10. plot excavation rate tunnel width vs. Gini
%*11. plot gini in vs gini out with multiple runs
%*12. plot cum. work vs. cum. pop
%*13. gini vs. gens for paper
%*14. gini vs reversal
%*15. ant exp lorenz with theory lorenz
%*16. plot q_bar vs rho_bar old type different ant num from single prob
%*17&18. V vs. N for 30 ant interpolated
%*19. code for another file to make movie from road data
%*20. excavation rate vs. time
%*21. plot actual ant exp densities
%*55. old gini vs generations
%************************************************************
showFigs=[21];
% showFigs=[7 16];
fold=uigetdir('D:\Projects\Ant_CA_GA\results');
filez=dir(fullfile(fold,'*.mat'));
NF=length(filez);
cc=jet(NF);

%% 1 gini vs gens for different numbers of ants
xx=1;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    
    
    fz= 20;
    % cc=get(gca,'colororder');
    cc=jet(NF);
    colormap(cc);
    
    
    gens=zeros(1,NF);
    ginis= cell(1,NF);
    ants=gens;
    maxGens=1;
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        gens(i)=gen;
        maxGens=max(maxGens,gen);
        ginis{i}=zeros(gen,2); %[generations, gini, ants]
        ginis{i}(:,1)=0:gen-1;
        ants(i)=nvars;
        for j=1:gen
            ginis{i}(j,2)=Gini(bestofgen{j});
        end
        plot(ginis{i}(:,1),ginis{i}(:,2),'-','linewidth',2,'color',cc(i,:))
        
    end
    [ants,inds]=sort(ants);
    ginis=ginis(inds);
    caxis([1 NF])
    
    cbarHandle = colorbar('YTick',...
        1+0.5*(NF-1)/NF:(NF-1)/NF:NF,...
        'YTickLabel',ants, 'YLim', [1 NF]); %[]=ants
    
    xlim([0,maxGens]);
    xlab=xlabel('Generation');ylab=ylabel('Gini Coefficient');
    title(cbarHandle,'N (ants)')
    figText(gcf,fz);
end
fold
%% 2 best of gen gini(last) vs. number of ants
xx=2;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    
    ginis=zeros(1,NF);
    ants=ginis;
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        ginis(i)=Gini(bestofgen{end});
        ants(i)=nvars;
    end
    [ants,inds]=sort(ants);
    ginis=ginis(inds);
    plot(ants,ginis,'o-','linewidth',2)
    
    % set(gca,'box','on','linewidth',2,'xtick',[0, 5,10, 15, 20], 'xticklabel',{'0','','10','','20'},'fontsize',fz);
    
    draw30line=1;
    if(draw30line)
        yval30=ginis(ants==30);
        plot([30,30],[0,yval30],'color','r','linewidth',2)
        plot([0,30],[yval30,yval30],'color','r','linewidth',2)
        set(gca,'box','on','linewidth',2,'xtick',[0, 25,30,50, 75, 100], 'xticklabel',{'0','','30','50','','100'});
        %         set(gca,'ytick',[0.2,0.4,0.6,yval30,0.8, 1.0], 'yticklabel',{'0.2','0.4','0.6',num2str(yval30,2),'0.8','1'});
    end
    % annotation('textarrow',[.225,.225],[.3,.8]);
    text(.75,.95,'unequal','fontsize',fz,'units','normalized');
    text(.75,.05,'equal','fontsize',fz,'units','normalized');
    xlabel('N (ants)');
    ylabel('Gini');
    figText(gcf,fz,2);
end
%% 3 number of ants, q vs rhos
xx=3;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    mark = 's';
    lw=2;
    fz= 20;
    [qm,rhom,ants]=deal(zeros(1,NF));
    cc=jet(NF);
    colormap(cc);
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        stopTime =3;
        startTime=.25;
        minz = 60*(stopTime-startTime); %%%uncomment
        %   minz = 60*(stopTime-startTime);
        stopTime = stopTime*7200;
        startTime=startTime*7200;
        ge=bestofgenOUT{end};
        prob=sort(ge/sum(ge));
        res=CA_FunctionsWill(prob,length(prob),432,TW,...
            energyMult,1,rechargeSteps,prob2turn,tuntip);
        singleLane = sum(res.occupied(startTime+1:stopTime,:),2)./TW; %/2 for 2 lanes
        
        if isfield(res,'density')
            dens = res.density(startTime+1:stopTime,:);
            d=1; %set zero for old way
        end
        singleFlow = sum(res.flow(startTime+1:stopTime,:),2)./TW;
        tunLen=(res.tunLength(startTime+1:stopTime,2));
        ts=120*minz; %timesteps to average over 120ts/min *60min/hour
        %     ss=size(res.occupied(1:stopTime,1));
        
        %         newq = reshape(singleFlow(1:end),ts,size(singleFlow(1:end),1)/ts);
        %         tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
        %         newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts);
        
        
        qq=sum(sum(res.markMatr(startTime:stopTime,2:end)))/res.pause2dig;
        newq=qq/ts;
        tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
        newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts)./tunLen;
        
        qm(i) = mean(newq);
        tunLen=mean(tunLen);
        if d
            rhom(i) = mean(dens);
        else
            rhom(i) = mean(newrho);
        end
        %
        %         rhom(i) = mean(dens);
        ants(i)=nvars;
    end
    %sort runs
    [ants,inds]=sort(ants);
    qm=qm(inds)*2;%timestep is .5 s
    rhom=rhom(inds);
    %plot newly sorted
    
    for(i=1:NF)
        plot(rhom(i),qm(i),'v','markerfacecolor',cc(i,:),'MarkerSize',15,'MarkerEdgeColor','k','LineWidth',lw);
    end
    
    haxis = gca;
    caxis([1 NF])
    aa=num2cell(ants);
    aa2=cellfun(@(x) num2str(x),aa,'uniformoutput',0);
    
    
    cbarHandle = colorbar('YTick',...
        1+0.5*(NF-1)/NF:(NF-1)/NF:NF,...
        'YTickLabel',{'5','','','20','','','','40','','','','60','','','','80','','','','100'}, 'YLim', [1 NF]); %[]=ants
    %     set(gca,'XTick',[0, 0.15 0.3]*1e-1,'YTick',[.25 .5 0.75]*1e-1);
    %     set(gca,'XTick',[0 0.005 0.01 0.015 0.02],'YTick',[.25 .5 0.75]*1e-1);
    
    %     set(cbarHandle,'fontsize',18');
    %     ylabel(cbarHandle,'N')
    title(cbarHandle,'N')
    %     axis([0 .4e-1 0, 8e-2])
    axis([0 .3e-1 0, 0.2])
    ax=gca;
    %     ax.YAxis.Exponent = -2;
    ax.XAxis.Exponent = -2;
    if(res.infEnergy)
        En=' E=\infty';
    else
        En=[];
    end
    
    %     set(gca,'box','on','linewidth',2);
    %     set(gcf,'position',[100,100,426,429]);
    xlab=xlabel('\rho (ants/BL) ');
    %     ylab=ylabel('q (ants/(BL\cdotmin)');
    ylab=ylabel('q (ants/s)');
    figText(gcf,fz,2);
end
%% 4 (gen(10)-gen(1)) vs. number of ants
xx=4;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    
    
    fz= 20;
    % cc=get(gca,'colororder');
    cc=jet(NF);
    colormap(cc);
    
    xlab=xlabel('Generation');ylab=ylabel('Gini Coefficient');
    gens=zeros(1,NF);
    ginis= zeros(NF,3);%[r1,rend,deltaG]
    ants=gens;
    maxGens=1;
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        ginis(i,1)=Gini(bestofgenOUT{1});
        ginis(i,2)=Gini(bestofgenOUT{end});
        ginis(i,3)=ginis(i,2)-ginis(i,1);
        ants(i)=nvars;
    end
    [ants,inds]=sort(ants);
    ginis=ginis(inds,:);
    plot(ants,abs(ginis(:,3)),'o-','linewidth',2)
    xlabel('number of ants');
    ylabel('\DeltaG');
    figText(gcf,fz,2);
end
%% *MAKES 5&6*
%Tunnel length excavated vs. time
%Tunnel length/ant excavated vs. time
xx=[5,6];
if(any(ismember(showFigs,xx)))
    figure(xx(1))
    hold on;
    fz= 20;
    
    for j=[0 1]%0=eq 1=uneq
        [V,ants]=deal(zeros(1,NF));
        for i=1:NF
            its=432;
            load(fullfile(fold,filez(i).name));
            if(j==0) %eq
                pp=ones(1,nvars);
                res=CA_FunctionsWill(pp,nvars,its,TW,...
                    energyMult,1,rechargeSteps,prob2turn,tuntip);
            else    %uneq
                %         res=CA_FunctionsWill(bestofgen{end},length(bestofgen{end}),numIts,TW,...
                res=CA_FunctionsWill(bestofgen{end},length(bestofgen{end}),432,TW,...
                    energyMult,1,rechargeSteps,prob2turn,tuntip);
            end
            
            pellet2grow=res.pell2grow;
            pellPerCm = pellet2grow*TW*2; %tunnwidth=2,another 2 for 1site=.5 cm
            pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerCm;
            ants(i)=nvars;
            V(i)=pellsTot/((its*100)/(3600*2));
        end
        %sort runs
        figure(xx(1));
        hold on;
        [ants,inds]=sort(ants);
        V=V(inds);
        %plot newly sorted
        plot(ants,V,'o-','linewidth',2);
        ylabel('V (cm/h)');
        xlabel('N (ants)');
        figText(gcf,fz,2);
        
        figure(xx(2));
        hold on;
        plot(ants,V./ants,'o-','linewidth',2);
        ylabel('V (cm/(h*ant)');
        xlabel('N (ants)');
        figText(gcf,fz,2);
    end
end

%% 7. *fixed* workload distribution q vs rho *use 16 for multires files
xx=7;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    mark = 's';
    lw=2;
    fz= 20;
    [qm,rhom,ants]=deal(zeros(1,NF));
    cc=jet(NF);
    colormap(cc);
    
    for i=1:NF %find ants =30
        load(fullfile(fold,filez(i).name));
        if nvars==30
            %%%%%is commented line diff?
            %             %     prob=bestofgen{end};
            prob= sort(bestofgen{end}/sum(bestofgen{end}));
        end
    end
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        stopTime =3;
        startTime=.25;
        minz = 60*(stopTime-startTime); %%%uncomment
        %   minz = 60*(stopTime-startTime);
        stopTime = stopTime*7200;
        startTime=startTime*7200;
        pp=InterpolateGAProbsFromProb(nvars,prob);
        
        %         res=CA_FunctionsWill(pp,length(pp),numIts,TW,...
        res=CA_FunctionsWill(pp,length(pp),432,TW,...
            energyMult,1,rechargeSteps,prob2turn,tuntip); %tuntip instead of 5
        singleLane = sum(res.occupied(startTime+1:stopTime,:),2)./TW; %/2 for 2 lanes
        
        if isfield(res,'density')
            dens = res.density(startTime+1:stopTime,:); %/2 for 2 lanes
            
        end
        singleFlow = sum(res.flow(startTime+1:stopTime,:),2)./TW;
        tunLen=(res.tunLength(startTime+1:stopTime,2));
        ts=120*minz; %timesteps to average over 120ts/min *60min/hour
        %     ss=size(res.occupied(1:stopTime,1));
        
        %         newq = reshape(singleFlow(1:end),ts,size(singleFlow(1:end),1)/ts);
        %         tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
        %         newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts);
        
        %get all pellets dug over start to stop time, number of pellets
        %represents 2x trips (1 to tunnel tip, and 1 back with pellet) we
        %can get successful flow rate this way
        qq=sum(sum(res.markMatr(startTime:stopTime,2:end)))/res.pause2dig;
        
        %(successful excavations)/(system TS *1s/2ts)= exc/s
        newq=qq/(ts/2);
        %each exc is 2 ant appearances: 2(ants/exc)*(exc/s)= (ants/s)
        newq=newq*2;
        tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
        
        newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts)./tunLen;%ants/(bL*lanes)
        
        
        qm(i) = mean(newq)/2; %%ants/(lanes*sec)
        tunLen=mean(tunLen);
        d=0; %set zero for old way
        if d
            rhom(i) = mean(dens);
        else
            rhom(i) = mean(newrho);
        end
        %
        %         rhom(i) = mean(dens);
        ants(i)=nvars;
    end
    %sort runs
    [ants,inds]=sort(ants);
    qm=qm(inds);
    rhom=rhom(inds);
    %plot newly sorted
    
    for(i=1:NF)
        plot(rhom(i),qm(i),'s','markerfacecolor',cc(i,:),'MarkerSize',15,'MarkerEdgeColor','k','LineWidth',lw);
    end
    
    haxis = gca;
    caxis([1 NF])
    aa=num2cell(ants);
    aa2=cellfun(@(x) num2str(x),aa,'uniformoutput',0);
    
    
    cbarHandle = colorbar('YTick',...
        1+0.5*(NF-1)/NF:(NF-1)/NF:NF,...
        'YTickLabel',{'5','','','20','','','','40','','','','60','','','','80','','','','100'}, 'YLim', [1 NF]); %[]=ants
    %     set(gca,'XTick',[0, 0.15 0.3]*1e-1,'YTick',[.25 .5 0.75]*1e-1);
    %     set(gca,'XTick',[0 0.005 0.01 0.015 0.02],'YTick',[.25 .5 0.75]*1e-1);
    
    %     set(cbarHandle,'fontsize',18');
    %     ylabel(cbarHandle,'N')
    title(cbarHandle,'N')
    %     axis([0 .4e-1 0, 8e-2])
    axis([0 .3e-1 0, 0.2])
    ax=gca;
    %     ax.YAxis.Exponent = -2;
    ax.XAxis.Exponent = -2;
    if(res.infEnergy)
        En=' E=\infty';
    else
        En=[];
    end
    
    %     set(gca,'box','on','linewidth',2);
    %     set(gcf,'position',[100,100,426,429]);
    xlab=xlabel('\rho (ants/(BL*lane)) ');
    %     ylab=ylabel('q (ants/(BL\cdotmin)');
    ylab=ylabel('q (ants/(s\cdotlane))');
    figText(gcf,fz,2);
end

%% 8 get cluster dissolution
xx=8;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    [V,ants]=deal(zeros(1,NF));
    h=waitbar(1/NF,[num2str(1),'/',num2str(NF)]);
    for i=1:NF
        tic
        if(i~=1)
            timeLeft=ttc*(NF-i);
            waitbar(i/NF,h,[num2str(i),'/',num2str(NF),' ',num2str(timeLeft,'%.0f'),'s left']);
        end
        load(fullfile(fold,filez(i).name));
        res=CA_FunctionsWill(bestofgen{end},length(bestofgen{end}),numIts,TW,...
            energyMult,1,rechargeSteps,prob2turn,tuntip);
        ants(i)=nvars;
        ttc=toc; %time to complete
    end
    closeWaitbar;
    %sort runs
    [ants,inds]=sort(ants);
    V=V(inds);
    %plot newly sorted
    plot(ants,V,'o-','linewidth',2);
    ylabel('');
    xlabel('');
    figText(gcf,fz,2);
end
%% 9 excavation amount for results for different single runs  V vs reversal
xx=9;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    cc=get(gca,'colororder');
    
    unitless=0;
    TW=2;
    datAll=zeros(NF,4);     %type,ant#,revProb,ExcavationAmount
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        [~,val]=parseFileNames(filez(i).name,1);
        xt=res.tunLength(:,1)*.5/3600;
        pellet2grow=res.pell2grow;
%         pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
%         pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerCm;
        pellPerBL = pellet2grow*TW; %tunnwidth=2,another 2 for 1site=.5 cm
        pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerBL;
        
        datAll(i,:)=[val(1),res.numants,val(2),pellsTot];
    end
    eq=datAll(datAll(:,1)==0,:); %equal distrib= type 0, square marker
    uneq=datAll(datAll(:,1)==1,:);%unequal distrib= type 1 circle marker
    
    % plot(eq(:,3),eq(:,4),'s-','markerfacecolor',cc(1,:),'MarkerSize',15,'MarkerEdgeColor','k','LineWidth',2);
    % plot(uneq(:,3),uneq(:,4),'o-','markerfacecolor',cc(2,:),'MarkerSize',15,'MarkerEdgeColor','k','LineWidth',2);
%     if(unitless)
%         TW=size(res.occupied,2);
%         ylabel('V/TW');
%         %         set(gca,'xtick',[0, .15 ,0.25 ,0.5,0.75 ,1], 'xticklabel',{'0','0.15','','0.5','','1'}, 'ytick',[0.5 1.0 1.5 1.82 2],'yticklabel', {'0.5','1','1.5','1.82','2'},'fontsize',fz);
%     else
%         TW=1;
%         ylabel('V (cm)');
%         %         set(gca,'xtick',[0, .15 ,0.25 ,0.5,0.75 ,1], 'xticklabel',{'0','0.15','','0.5','','1'}, 'ytick',[1 2 3 3.642 4],'yticklabel', {'1','2','3','3.64','4'},'fontsize',fz);
%     end
%     
    %%%%%with markers%%%
    %     plot(eq(:,3),eq(:,4)/TW,'s-','markerfacecolor',cc(1,:),'MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1.4);
    %     plot(uneq(:,3),uneq(:,4)/TW,'o-','markerfacecolor',cc(2,:),'MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1.4);
    
    %%%%%without markers%%%%%%%%
    plot(eq(:,3),eq(:,4),'-','LineWidth',3);
    plot(uneq(:,3),uneq(:,4),'-','LineWidth',3);
    
    [mRev,revInd]=max(uneq(:,4));
%     pts('mRev=',mRev,' reversal=',uneq(revInd,3));
    
    plot([uneq(revInd,3),uneq(revInd,3)],[0,mRev]/TW,'r','linewidth',2);
    plot([0,uneq(revInd,3)],[mRev,mRev]/TW,'r','linewidth',2);
    
    
    
    xlabel('Reversal probability');
    ylabel('V/TW (BL/BW)');
    figText(gcf,fz,2);
    
    set(gca,'xtick',[0 .1 .25 .5 .75 1])
    set(gca,'ytick',[ 5 10 15 20 25 30],'yticklabel',{'','10','','20','','30'})
    % set(gca,'yscale','log','xscale','log','xticklabelmode','auto','xtickmode','auto','yticklabelmode','auto','ytickmode','auto')
end

%% 10. plot excavation rate tunnel width vs. Gini
%to get 3 different tunnel widths edit TW to be 2 3 4 etc
xx=10;
if(showFigs(showFigs==xx))
    figure(xx(1))
    hold on;
    fz= 20;
    for qq=[2,3,4]
    TW=qq;
    load(fullfile(pwd,'gini.mat'));
    gDatAll=zeros(length(giniX),4);
    tt=432*4;
    parfor i=1:length(giniX)
        res=CA_FunctionsWill(giniY(i,:),length(giniY(i,:)),tt,TW,...
            0,0,600,.3,10);
        pellet2grow=res.pell2grow;
%         pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
        pellPerBL = pellet2grow*TW; %tunnwidth=2,took out extra 2 because we are plotting vs. BL
        pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerBL;
        Gout=Gini(sum(res.markMatr(:,2:end)));
        gDatAll(i,:)=[res.numants,pellsTot,G(i),Gout];
    end
    %with time in denom
%     plot(gDatAll(:,3),gDatAll(:,2)*3600/((tt/2)*100),'o-','markerfacecolor','w','linewidth',2,'markersize',7);
%     without time in denom 
    plot(gDatAll(:,3),gDatAll(:,2),'o-','markerfacecolor','w','linewidth',2,'markersize',7);
        xlabel('Gini');
%     ylabel('V/TW (BL/(BW*h))');
        ylabel('V/TW (BL/BW)');
    %     set(gca,'xtick',[0, .15 ,0.25 ,0.5,0.75 ,1], 'xticklabel',{'0','0.15','','0.5','','1'}, 'ytick',[1 2 3 3.642 4],'yticklabel', {'1','2','3','3.64','4'},'fontsize',fz);
    %     set(gca,'yscale','log','xscale','log','xticklabelmode','auto','xtickmode','auto','yticklabelmode','auto','ytickmode','auto')
    
    figText(gcf,fz,2);
    xlim([0, 1]);
    yl=ylim;
    ylim([yl(1)*.975,yl(2)*1.025]);
    set(gca,'xtick',[0, .25 ,0.5 ,.75 ,1]);
    end
end

%% 11. plot gini in vs gini out with multiple runs
xx=11;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    numIts=5;
    TW=2;
    load(fullfile(pwd,'gini.mat'));
    Gout=zeros(length(giniX),numIts);
    for i=1:length(giniX)
        for j=1:numIts
            res=CA_FunctionsWill(giniY(i,:),length(giniY(i,:)),432,TW,...
                1,1,10,.3,10);

            pellet2grow=res.pell2grow;
            pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
            pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerCm;
            Gout(i,j)=Gini(sum(res.markMatr(:,2:end)));
        end
        pts(i,'/',length(giniX));
    end
    xlabel('Gini');
    ylabel('V (cm)');
    %     set(gca,'xtick',[0, .15 ,0.25 ,0.5,0.75 ,1], 'xticklabel',{'0','0.15','','0.5','','1'}, 'ytick',[1 2 3 3.642 4],'yticklabel', {'1','2','3','3.64','4'},'fontsize',fz);
    %     set(gca,'yscale','log','xscale','log','xticklabelmode','auto','xtickmode','auto','yticklabelmode','auto','ytickmode','auto')
    mGout=mean(Gout,2);errGout=std(Gout,0,2);
    errorbar(giniX,mGout,errGout,'linewidth',2);
    plot([0,1],[0,1],':k','linewidth',3);
    xlabel('G_{in}');
    ylabel('G_{out}');
    figText(gcf,fz,2);
    xlim([0,1]);
    yl=ylim;
    ylim([yl(1)*.975,yl(2)*1.025]);
    set(gca,'xtick',[0, .2 ,.4 ,.6 ,.8 ,1],'ytick',[0, .2 ,.4 ,.6 ,.8 ,1]);
    axis square;
end
%% 12. plot cum. work vs. cum. pop
xx=12;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    %     load(fullfile(pwd,'results','v1.mat'));
    %     y=bestofgen{end};
    TW=2;
    co =[184 77 157]./255;
    textz = 'mode 3';
    %uncomment if not using matrix file containing correct tunnel tip=3
    %     res =CA_Functions2(y,length(y),432,2,1,1,600,.3,3);  %probs,numants,numits*10000,width,infEnergy
    %     load(['D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge\long rand 50 gens\N=30_tw=2_2017-10-16-02-45.mat']);
    load(fullfile(fold,filez(1).name));
    GENNUM=length(bestofgen);%generation number to plot from
    y=bestofgen{GENNUM};
    res=CA_FunctionsWill(y,length(y),numIts,TW,energyMult,0,rechargeSteps,prob2turn,tuntip);
    %     res=CA_FunctionsWill(y,length(y),432,TW,1,1,10,.3,10);
    z = zeros(1,29);
    z(30)=1;
    [~,eqLine]=Gini(ones(1,30));
    [~,uneqLine]=Gini(z);
    [~,theoryLine]=Gini(sum(res.markMatr(:,2:end)));
    cc=get(gca,'colororder');
    ucol=cc(1,:);
    ecol=cc(2,:);
    tcol=cc(5,:);
    
    plot(eqLine(:,1),eqLine(:,2),'linewidth',3,'color',ecol);
    plot(uneqLine(:,1),uneqLine(:,2),'linewidth',3,'color',ucol);
    plot(theoryLine(:,1),theoryLine(:,2),'linewidth',3,'color',tcol);
    fz=20;
    
    text(.35,.61,'Equal','fontsize',fz,'color',ecol,'fontname','arial','fontangle','italic');
    text(.45,.325,'Theory','fontsize',fz,'color',tcol,'fontname','arial','fontangle','italic');
    text(.55,.0725,'Unequal','fontsize',fz,'color',ucol,'fontname','arial','fontangle','italic');
    
    xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}');
    set(gca,'box','on','linewidth',2,'xtick',[0 .5 1], 'ytick',[.5 1],'fontsize',fz);
    axis([0 1 0 1]);
    axis square
end
%% 13 gini vs generations for paper
xx=13;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    cc=get(gca,'colororder');
    fz= 20;
    gens=0:50;
    %initialize gini vars
    [giniEq,giniUneq,giniRand]=deal(gens);
    
    
    %     %unequal seed
    uneq=dir(fullfile(fold,'uneq','*.mat'));
    load(fullfile(uneq.folder,uneq.name));
    for i = 1:length(gens)
        [giniUneq(i),~]=Gini(bestofgen{i});%out goes down
    end
    
    
    %equal seed
    eq=dir(fullfile(fold,'eq','*.mat'));
    load(fullfile(eq.folder,eq.name));
    for i = 1:length(gens)
        [giniEq(i),~]=Gini(bestofgen{i}); %out goes down
        %         res=CA_FunctionsWill(bestofgen{i},length(bestofgen{i}),432,TW,1,1,10,.3,10);
        %         [giniEq(i),~]=Gini(sum(res.markMatr(:,2:end)));
    end
    
    %rand seed
    ra=dir(fullfile(fold,'rand','*.mat'));
    load(fullfile(ra.folder,ra.name));
    for i = 1:length(gens)
        [giniRand(i),~]=Gini(bestofgen{i});%out goes down
    end
    %
    plot(gens,giniEq   ,'-','linewidth',3,'color',cc(1,:));
    plot(gens,giniUneq ,'-','linewidth',3,'color',cc(2,:));
    plot(gens,giniRand ,'-','linewidth',3,'color',cc(5,:));
    
    axis([0 length(gens)-1 0 1]);
    
    set(gca,'box','on','linewidth',2);
    set(gca,'xtick',[0,25, 50], 'xticklabel',{'0','','50'} ,'fontsize',fz);
    %     'ytick',[ .25,.31, .5 .75 1],'yticklabel', {'0.31','0.5','','1'}
    axis square
    
    xlab=xlabel('Generation');ylab=ylabel('Gini Coefficient');
    figText(gcf,fz);
end


%% 14 gini vs reversal
xx=14;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    cc=get(gca,'colororder');
    
    unitless=1;
    
    datAll=zeros(NF,4);     %type,ant#,revProb,ExcavationAmount
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        [~,val]=parseFileNames(filez(i).name,1);
        a=Gini(sum(res.markMatr(:,2:end)));
        datAll(i,:)=[val(1),res.numants,val(2),a];
    end
    eqInds=find(datAll(:,1)==0);
    uneqInds=find(datAll(:,1)==1);
    plot(datAll(uneqInds,3),datAll(uneqInds,4),'o-','linewidth',2);
    plot(datAll(eqInds,3),datAll(eqInds,4),'o-','linewidth',2);
    
    
    
    xlabel('Reversal probability');
    ylabel('Gini');
    figText(gcf,fz,2);
    % set(gca,'yscale','log','xscale','log','xticklabelmode','auto','xtickmode','auto','yticklabelmode','auto','ytickmode','auto')
end
%% 15. ant exp lorenz with theory lorenz
xx=15;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    clear xy
    TW=2;
    % load([pwd,'\R=600P=0.45\GAdat.mat']);
    %     load('D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge\long rand 50 gens\N=30_tw=2_2017-10-16-02-45.mat');
    load(fullfile(fold,filez(1).name));
    runs = 2;
    GEN=length(bestofgen);
    y=bestofgen{GEN};
    % % %     for i=1:runs
    % % %         %prob,ants,time,tw,energymult,0,rech,ptt,tuntip
    % % %
    % % %         res=CA_FunctionsWill(y,length(y),numIts,TW,energyMult,0,rechargeSteps,prob2turn,tuntip);
    % % %         [ggruns(i),xy{i}]=Gini(sum(res.markMatr(:,2:end)));
    % % %         ggruns
    % % % %                 [ggruns(i),xy{i}]=Gini(bestofgen{GEN});
    % % %         %              CA_FunctionsWill(prob,length(prob),ni,2,0,1,rec,ptt);
    % % %     end
    % % %     yy=cell2mat(cellfun(@(x) x(:,2),xy,'UniformOutput',0));
    % % %     meanY=mean(yy,2);
    % % %     % meanY=mean([a{:}],2);
    % % %     meanG=mean(ggruns);
    % % %     err=std(yy,1,2);
    % % %     cc=parula(15);
    % % %     shadedErrorBar(xy{1}(:,1),meanY,err,{'Color',cc(10,:),'LineWidth',3},1);
    [gx,gy]=Gini(y);
    
    % text(.651,.829,['Theory'],'Color',cc(10,:),'fontsize',fz);
    % xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}');
    % figText(gcf,18);
    % xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}');
    set(gca,'box','on','linewidth',2);
    set(gcf,'Position',[1 1 594 562]);
    
    
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
    axis([0 1 0 1]);
    
    plot(gy(:,1),gy(:,2),'Linewidth',3,'color',cc(9,:));
    figText(gcf,fz);
end

%% 16. plot q_bar vs rho_bar old type different ant num from single prob
xx=16;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    titl='Traffic (CA)';
    %     load(fullfile(fold,filez.name)); %load up so I can get bestofgens var
    % eq=0; %% zero for inequal work
    fz=20;
    TW=2;
    for eq = [0 1]
        
        if(eq)
            %     fold ='\eqResSavesInfEnergy\';
            %     fold ='\eqResSavesInfEnergyPell600\';
            folddz=fullfile(fold,'diffNeq');
            %             fold = ['D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\diffNeq'];
            mark = 's';
            lw=2;
        else
            %     fold ='\uneqResSavesInfEnergy\';
            %     fold ='\uneqResSavesInfEnergyPell600\';
            %             fold =['D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\diffNuneq'];
            folddz=fullfile(fold,'diffNuneq');
            mark = 'o';
            lw = 1.9;
        end
        
        
        flist=dir([folddz,'\*.mat']);
        for i = 1:length(flist)
            num(i)=sscanf(flist(i).name,'res%d');
        end
        
        
        
        
        
        
        [~,order] = sort(num);
        flist(:)=flist(order);
        ants=[2 3 5 8 10 12 15 18 20 25 30 50 100];
        numcolors = length(ants*3);
        cmap = jet(numcolors);
        stopTime =3;
        startTime=.25;
        minz = 60*(stopTime-startTime); %%%uncomment
        %   minz = 60*(stopTime-startTime);
        stopTime = stopTime*7200;
        startTime=startTime*7200;
        for i =1:length(flist)
            flist(i).name;
            load(fullfile(folddz,flist(i).name));
            stopTime =3;
            startTime=.25;
            minz = 60*(stopTime-startTime); %%%uncomment
            %   minz = 60*(stopTime-startTime);
            stopTime = stopTime*7200;
            startTime=startTime*7200;
            %         pp=InterpolateGAProbsFromProb(nvars,prob);
            %
            % %         res=CA_FunctionsWill(pp,length(pp),numIts,TW,...
            %         res=CA_FunctionsWill(pp,length(pp),432,TW,...
            %             energyMult,1,rechargeSteps,prob2turn,tuntip); %tuntip instead of 5
            
            singleLane = sum(res.occupied(startTime+1:stopTime,:),2)./TW; %/2 for 2 lanes
            
            if isfield(res,'density')
                dens = res.density(startTime+1:stopTime,:); %/2 for 2 lanes
                
            end
            singleFlow = sum(res.flow(startTime+1:stopTime,:),2)./TW;
            tunLen=(res.tunLength(startTime+1:stopTime,2));
            ts=120*minz; %timesteps to average over 120ts/min *60min/hour
            %     ss=size(res.occupied(1:stopTime,1));
            
            %         newq = reshape(singleFlow(1:end),ts,size(singleFlow(1:end),1)/ts);
            %         tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
            %         newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts);
            
            %get all pellets dug over start to stop time, number of pellets
            %represents 2x trips (1 to tunnel tip, and 1 back with pellet) we
            %can get successful flow rate this way
            qq=sum(sum(res.markMatr(startTime:stopTime,2:end)))/res.pause2dig;
            %%%%%%%%%%%%%
            newq=qq/(ts/2); %each exc is 2 ant appearances: 2(ants/exc)*(exc/s)= (ants/s)     
            qm(i) = newq/TW; %%ants/(lanes*sec)
            
            %%%%%%%%%%%%%%
            tunLen=reshape(tunLen(1:end),ts,size(tunLen(1:end),1)/ts);
            
%             newrho=reshape(singleLane(1:end),ts,size(singleLane(1:end),1)/ts)./tunLen;%ants/(bL*lanes)
            newrho=mean(singleLane./tunLen);
            
            
%             tunLen=mean(tunLen);
            d=1; %set zero for old way
            if d
                rhom(i) = mean(dens.*tunLen);
            else
                rhom(i) = mean(newrho);
            end
            %
            %         rhom(i) = mean(dens);
            %         ants(i)=nvars;
            ants(i)=res.numants;
        end
        
        %sort runs
        [ants,inds]=sort(ants,'ascend');
        qm=qm(inds);
        rhom=rhom(inds);
        %plot newly sorted
        
        for(i=1:length(flist))
            plot(rhom(i),qm(i),mark,'markerfacecolor',cmap(i,:),'MarkerSize',15,'MarkerEdgeColor','k','LineWidth',lw);
        end
        
        
        haxis = gca;
        set(gca,'fontsize',55);
        figText(gcf,16)
        colormap(cmap);
        
        caxis([1 numcolors])
        
        cbarHandle = colorbar('YTick',...
            1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors,...
            'YTickLabel',ants, 'YLim', [1 numcolors]); %[]=ants
        
        %         set(gca,'XTick',[0.25 0.5],'YTick',[0 .1 .2],'XLim',[0.045 0.55],'YLim',[0,0.15]);
        %             set(gca,'XTick',[0.025 0.05 0.075],'YTick',[0.05 0.1 0.15],'yTicklabel',[], 'xTicklabel',[],'fontsize',fz);
%         set(gca,'XTick',[0.25 0.5 0.75],'YTick',[0.05 0.1 0.15],'fontsize',fz);
        
        %     set(cbarHandle,'fontsize',18');
        %     axis([0 .4 0 .075])
        if(res.infEnergy)
            En=' E=\infty';
        else
            En=[];
        end
        % title([titl,' (',num2str(stopTime/7200),'h)',En]);
        % title([titl,En]);
        set(gca,'box','on','linewidth',2);
        set(gca,'position',[ 0.1380    0.1173    0.6301    0.8077]);
        set(gcf,'position',[   680   495   576   483]);
        xlab=xlabel('\rho (ants/(BL*lane)) ');
        %     ylab=ylabel('q (ants/(BL\cdotmin)');
        ylab=ylabel('q (ants/(s\cdotlane))');
%         ylim([0,.15]);
        xlim([0,.9])
        ylim([0,.1])
        set(gca,'ytick',[.02 .04 .06 .08 0.1])
    end
    
    a = 0.001:.001:.3; % alpha
g = 1/25; % "switching rate" of excavation 
v = .8; % Speed of ant (1 cell per frame)
L = 5; % Length of tunnel
S = .8; % Switching rate due to reversal
p = 2*((L/2-.5*(1./a-v/g))./(.5*(1./a-v/S))+2)/v.*(1./(v*a)+2*L/v+1/S.*((L-.5.*(1./a-v/g))./(.5*(1./a-v/S))+1)+1/g).^-1;
% p is rho, or density of tunnel

q = (1./(v*a)+2*L/v+1/S.*((L-.5.*(1./a-v/g))./(.5*(1./a-v/S))+1)+1/g).^-1;
% q is flow rate, in ants per frame

conv = 0.5; % conversion, 0.5 seconds per frame.

plot(p,2*q/conv,'linewidth',4);

end
%% *MAKES 17&18* V vs. N for 30 ant interpolated
%Tunnel length excavated vs. time
%Tunnel length/ant excavated vs. time
xx=[17,18];
if(any(ismember(showFigs,xx)))
    figure(xx(1))
    hold on;
    fz= 20;
    ants = [2 3 5 8 10 12 15 18 20 25 30 50 100];
    for j=[0 1]%0=eq 1=uneq
        V=zeros(1,length(ants));
        for i=1:length(ants)
            load(['D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\1-100_n_ants\N=30_tw=2_2017-10-22-23-01.mat']);
            if(j==0) %eq
                pp=ones(1,ants(i));
                res=CA_FunctionsWill(pp,ants(i),432,TW,...
                    energyMult,1,rechargeSteps,prob2turn,tuntip);
            else    %uneq
                pp=InterpolateGAProbsFromProb(ants(i),bestofgen{end});
                %         res=CA_FunctionsWill(bestofgen{end},length(bestofgen{end}),numIts,TW,...
                res=CA_FunctionsWill(pp,length(pp),432,TW,...
                    energyMult,1,rechargeSteps,prob2turn,tuntip);
            end
            xt=res.tunLength(:,1)*.5/3600;
            pellet2grow=res.pell2grow;
            pellPerCm = pellet2grow*TW*2; %tunnwidth=2,another 2 for 1site=.5 cm
            pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerCm;
            %             ants(i)=nvars;
            V(i)=pellsTot/(43200/(3600*2));
        end
        %sort runs
        figure(xx(1));
        hold on;
        [ants,inds]=sort(ants);
        V=V(inds);
        %plot newly sorted
        plot(ants,V,'o-','linewidth',2);
        ylabel('V (cm/h)');
        xlabel('N (ants)');
        figText(gcf,fz,2);
        
        figure(xx(2));
        hold on;
        plot(ants,V./ants,'o-','linewidth',2);
        ylabel('V (cm/(h*ant)');
        xlabel('N (ants)');
        figText(gcf,fz,2);
    end
end

%% 19 code for another file to make movie from road data
xx=19;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    set(gcf,'position',[-177,735,1957,84])
    load(fullfile(pwd,'gini.mat'));
    [res,rf]=CA_FunctionsWill(giniY(8,:),length(giniY(8,:)),432,TW,...
        0,0,600,.3,10);
    rf2=rf+4;
    colormap(lines)
    c=colormap;
    c(4,:)=[0,0,0];
    c(5,:)=[1,1,1];
    colormap(c);
    for i=100:300
        clf;
        image(rf2(:,:,i))
        pause(0.5);
    end
end
%% 20 excavation rate vs. time
xx=20;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    TW=2;
%     set(gcf,'position',[-177,735,1957,84])
    load('D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\1-100_n_ants\N=30_tw=2_2017-10-22-23-01.mat');
    tt=432*4*100;
    prob=bestofgen{end};
    
    res=CA_FunctionsWill(prob,length(prob),tt/100,TW,...
        0,0,600,.3,10);
    pellet2grow=res.pell2grow;
%         pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
    pellPerBL = pellet2grow*TW; %tunnwidth=2,took out extra 2 because we are plotting vs. BL
    pellsTot = cumsum(sum(res.markMatr(:,2:end),2)/res.pause2dig/pellPerBL);
    pellsTot = pellsTot(1:7200:end);
    
    simTime=[0:length(pellsTot)-1];
    plot(simTime,pellsTot,'-o','markerfacecolor','w','linewidth',2)
    xlabel('time (h)');
    ylabel('BL');
    figText(gcf,16); 
    
end
%% 21. plot actual ant exp densities
xx=21;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    xl=[0,4]; %xlim
    yl=[0,1]; %ylim
%%%%%%%%%%%%%%%old data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     load('antno.mat');
%     antBL=6; %6mm/BL
%     antBW=2;
%     TLants=TL/antBL;% tunnel length/ant length gives TL in units of BL
%     e1=mean([exp1120(1,:),exp1120(2,:),exp1120(3,:)])./TLants(1)/antBW;
%     e2=mean([exp1123(1,:),exp1123(2,:),exp1123(3,:)])./TLants(2)/antBW;
%     e3=mean([exp1130(1,:),exp1130(2,:),exp1130(3,:)])./TLants(3)/antBW;
%     plot([1 2 3],[e1 e2 e3],'-o');   
%     
%     meanDat=mean([e1,e2,e3]);
%     plot(xl,[meanDat,meanDat],'r','linewidth',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load('antnoNew.mat');
    antBL=6; %6mm/1BL
    antBW=2; %2BW/1 tunnel
    
    inds = [4:6];
    TL1ants=antBW*TL(1,:)./(antBL); % (BW/tunnel)*mm/(mm/BL)=BW*BL/tunn
    TL2ants=antBW*TL(2,:)./(antBL);
    TL3ants=antBW*TL(3,:)./(antBL);
    
    e1=exp1120./TL1ants'; %(ants/tunnel)/(BL*BW/tunnel)= ants/(BL*BW)
    e2=exp1123./TL2ants';
    e3=exp1130./TL3ants';
    
    e1=e1(inds,:);
    e2=e2(inds,:);
    e3=e3(inds,:);
    
    e1m=mean(mean(e1,2));
    e2m=mean(mean(e2,2));
    e3m=mean(mean(e3,2));
    
    e1Err=std(mean(e1,2));
    e2Err=std(mean(e2,2));
    e3Err=std(mean(e3,2));
    
    errorbar([1 2 3],[e1m e2m e3m],[e1Err e2Err, e3Err]);
    meanDat=mean([e1m,e2m,e3m]);
    plot(xl,[meanDat,meanDat],'-','linewidth',2);

    
    set(gca,'xtick',[1 2 3], 'xticklabels',{'EXP 1', 'EXP 2', 'EXP 3'})
%     set(gca,'ytick',[.2  .4 meanDat .6 .8 1], 'yticklabels',{'0.2','0.4',num2str(meanDat,2),'0.6','0.8','1' })
    
    xlim(xl);
    ylim(yl);
    xlabel('experiment');
    ylabel('\rho (ants/(BL\cdotBW)');
    figText(gcf,16);
    
    
end
%% 55 old gini vs generations
xx=55;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    
    
    %equal seed
    gens20=0:20;
    load(['D:\Projects\Ant_CA_GA\results\largerTunnWidth\100antstunWidth=5.mat']);
    load(['D:\Projects\Ant_CA_GA\results\largerTunnWidth\100antstunWidth=5.mat']);
    % gens=0:length(bestofgenOUT)-1;
    
    gini100=zeros(1,length(gens20));
    for i = 1:length(gens20)
        [gini100(i),~]=Gini(bestofgenOUT{i});
    end
    
    
    %unequal seed
    gens10=0:10;
    load('D:\Projects\Ant_CA_GA\results\largerTunnWidth\150antstunWidth=5.mat');
    gini150=zeros(1,length(gens10));
    for i = 1:length(gens10)
        [gini150(i),~]=Gini(bestofgenOUT{i});
    end
    
    
    
    
    plot(gens20,gini100     ,'-','linewidth',3,'color',cc(1,:));
    plot(gens10,gini150   ,'-','linewidth',3,'color',cc(2,:));
    % plot(gens,randGini   ,'-','linewidth',3,'color',cc(5,:));
    text(9,.2,'100 ants'  ,'color',cc(1,:),'fontsize',fz);
    text(9,.8,'150 ants','color',cc(2,:),'fontsize',fz);
    % text(9,.4,'Random' ,'color',cc(5,:),'fontsize',fz);
    axis([0 length(gens20)-1 0 1]);
    % set(gca,
    % set(gcf,'resize','off');
    set(gca,'box','on','linewidth',2,'xtick',[0, 5,10, 15, 20], 'xticklabel',{'0','','10','','20'}, 'ytick',[ 0,.25,.31, .5 .75 1],'yticklabel', {'0','','0.31','0.5','','1'},'fontsize',fz);
    % set(xlab,'units','normalized','Position',[.5 -.02 0]);
    % set(ylab,'units','normalized','Position',[0 .5 0]);
    % set(xlab,'position',get(xlab,'position')+[25 .2  0],'margin',1);
    axis square
    % set(xlab,'margin',3,'clipping','on','position',get(xlab,'position')+[25 .01  0]);
end
