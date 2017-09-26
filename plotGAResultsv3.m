clear all
% close all
% fold=uigetdir('A:\SmarticleRun\')
% load(fullfile(fold,'amoebaData.mat'));

% load('A:\SmarticleRun\Amoeba_newsquare_1_dead\amoebaData.mat');



%************************************************************
%* Fig numbers:
%* 1. gini vs gens
%* 2. best of gen gini(last) vs. number of ants
%* 3. number of ants, q vs rho *use 7*
%* 4. (gen(10)-gen(1)) vs. number of ants
%* 5/6. Tunnel length excavated vs. time  *AND* per ant
%* 7. *fixed* workload distribution for q vs rho
%* 8. get cluster dissolution
%* 9. excavation amount for results for different single runs
%*10. plot excavation rate vs. Gini in
%*11. plot gini in vs gini out with multiple runs
%*15. template
%************************************************************
showFigs=[10];
fold=uigetdir('D:\Projects\Ant_CA_GA\results');
filez=dir(fullfile(fold,'*.mat'));
NF=length(filez);
cc=jet(NF);
%% 1 gini vs gens
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
        ginis(i)=Gini(bestofgenOUT{end});
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
        set(gca,'ytick',[0.2,0.4,0.6,yval30,0.8, 1.0], 'yticklabel',{'0.2','0.4','0.6',num2str(yval30,2),'0.8','1'});
    end
    % annotation('textarrow',[.225,.225],[.3,.8]);
    text(.75,.95,'unequal','fontsize',fz,'units','normalized');
    text(.75,.05,'equal','fontsize',fz,'units','normalized');
    xlabel('number of ants');
    ylabel('Gini coefficient');
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
        prob=sort(bestofgenOUT{end}/sum(bestofgenOUT{end}));
        res=CA_FunctionsWill(prob,length(prob),numIts,TW,...
            energyMult,1,rechargeSteps,prob2turn,tuntip);
        singleLane = sum(res.occupied(startTime+1:stopTime,:),2)./TW; %/2 for 2 lanes
        
        if isfield(res,'density')
            dens = res.density(startTime+1:stopTime,:); %/2 for 2 lanes
            d=0; %set zero for old way
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
    [V,ants]=deal(zeros(1,NF));
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        res=CA_FunctionsWill(bestofgen{end},length(bestofgen{end}),numIts,TW,...
            energyMult,1,rechargeSteps,prob2turn,tuntip);
        xt=res.tunLength(:,1)*.5/3600;
        pellet2grow=res.pell2grow;
        pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
        pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerCm;
        ants(i)=nvars;
        V(i)=pellsTot;
    end
    %sort runs
    [ants,inds]=sort(ants);
    V=V(inds);
    %plot newly sorted
    plot(ants,V,'o-','linewidth',2);
    ylabel('V (cm)');
    xlabel('N (ants)');
    figText(gcf,fz,2);
    
    figure(xx(2));
    hold on;
    plot(ants,V./ants,'o-','linewidth',2);
    ylabel('V (cm/ant)');
    xlabel('N (ants)');
    figText(gcf,fz,2);
end

%% 7 fixed workload distribution for q vs rho
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
            %     prob=bestofgen{end};
            prob= sort(bestofgenOUT{end}/sum(bestofgenOUT{end}));
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
        res=CA_FunctionsWill(pp,length(pp),numIts,TW,...
            energyMult,1,rechargeSteps,prob2turn,tuntip); %tuntip instead of 5
        singleLane = sum(res.occupied(startTime+1:stopTime,:),2)./TW; %/2 for 2 lanes
        
        if isfield(res,'density')
            dens = res.density(startTime+1:stopTime,:); %/2 for 2 lanes
            d=0; %set zero for old way
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
        rhom(i) = mean(dens);
        ants(i)=nvars;
    end
    %sort runs
    [ants,inds]=sort(ants);
    qm=qm(inds)*2;%timestep is .5 s
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
    xlab=xlabel('\rho (ants/BL) ');
    %     ylab=ylabel('q (ants/(BL\cdotmin)');
    ylab=ylabel('q (ants/s)');
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
%% 9 excavation amount for results for different single runs
xx=9;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    cc=get(gca,'colororder');
    datAll=zeros(NF,4);     %type,ant#,revProb,ExcavationAmount
    for i=1:NF
        load(fullfile(fold,filez(i).name));
        [~,val]=parseFileNames(filez(i).name,1);
        xt=res.tunLength(:,1)*.5/3600;
        pellet2grow=res.pell2grow;
        pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
        pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerCm;
        datAll(i,:)=[val(1),res.numants,val(2),pellsTot];
    end
    eq=datAll(datAll(:,1)==0,:); %equal distrib= type 0, square marker
    uneq=datAll(datAll(:,1)==1,:);%unequal distrib= type 1 circle marker
    
    % plot(eq(:,3),eq(:,4),'s-','markerfacecolor',cc(1,:),'MarkerSize',15,'MarkerEdgeColor','k','LineWidth',2);
    % plot(uneq(:,3),uneq(:,4),'o-','markerfacecolor',cc(2,:),'MarkerSize',15,'MarkerEdgeColor','k','LineWidth',2);
    
    plot(eq(:,3),eq(:,4),'s-','markerfacecolor',cc(1,:),'MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1.4);
    plot(uneq(:,3),uneq(:,4),'o-','markerfacecolor',cc(2,:),'MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1.4);
    
    [mRev,revInd]=max(uneq(:,4));
    pts('mRev=',mRev,' reversal=',uneq(revInd,3));
    
    plot([uneq(revInd,3),uneq(revInd,3)],[0,mRev],'r','linewidth',2);
    plot([0,uneq(revInd,3)],[mRev,mRev],'r','linewidth',2);
    
    set(gca,'xtick',[0, .15 ,0.25 ,0.5,0.75 ,1], 'xticklabel',{'0','0.15','','0.5','','1'}, 'ytick',[1 2 3 3.642 4],'yticklabel', {'1','2','3','3.64','4'},'fontsize',fz);
    ylabel('V (cm)');
    xlabel('Reversal probability');
    figText(gcf,fz,2);
    % set(gca,'yscale','log','xscale','log','xticklabelmode','auto','xtickmode','auto','yticklabelmode','auto','ytickmode','auto')
end

%% 10 plot excavation rate vs. Gini
xx=[10];
if(showFigs(showFigs==xx))
    figure(xx(1))
    hold on;
    fz= 20;
    TW=4;
    load(fullfile(pwd,'gini.mat'));
    gDatAll=zeros(length(giniX),4);
    for i=1:length(giniX)
        res=CA_FunctionsWill(giniY(i,:),length(giniY(i,:)),432,TW,...
            1,1,10,.3,10);
        xt=res.tunLength(:,1)*.5/3600;
        pellet2grow=res.pell2grow;
        pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
        pellsTot = sum(sum(res.markMatr(:,2:end)))/res.pause2dig/pellPerCm;
        Gout=Gini(sum(res.markMatr(:,2:end)));
        gDatAll(i,:)=[res.numants,pellsTot,G(i),Gout];
    end
    %Gini In
    plot(gDatAll(:,3),gDatAll(:,2)/TW,'o-','markerfacecolor','w','linewidth',2,'markersize',7);
    %Gini Out
%     plot(gDatAll(:,4),gDatAll(:,2),'o-','markerfacecolor','w','linewidth',2,'markersize',7);
    xlabel('Gini');
    ylabel('V/TW');
    %     set(gca,'xtick',[0, .15 ,0.25 ,0.5,0.75 ,1], 'xticklabel',{'0','0.15','','0.5','','1'}, 'ytick',[1 2 3 3.642 4],'yticklabel', {'1','2','3','3.64','4'},'fontsize',fz);
    %     set(gca,'yscale','log','xscale','log','xticklabelmode','auto','xtickmode','auto','yticklabelmode','auto','ytickmode','auto')
    
    figText(gcf,fz,2);
    xlim([0, 1]);
    yl=ylim;
    ylim([yl(1)*.975,yl(2)*1.025]);
    set(gca,'xtick',[0, .2 ,.4 ,.6 ,.8 ,1]);

end

%% 11. plot gini in vs gini out with multiple runs
xx=[11];
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fz= 20;
    numIts=5;
    TW=4;
    load(fullfile(pwd,'gini.mat'));
    Gout=zeros(length(giniX),numIts);
    for i=1:length(giniX)
        parfor j=1:numIts
        res=CA_FunctionsWill(giniY(i,:),length(giniY(i,:)),432,TW,...
            1,1,10,.3,10);
        xt=res.tunLength(:,1)*.5/3600;
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
