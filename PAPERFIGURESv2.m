%% OPTIMAL LORENZ(*PAPER QUALITY SET TUNNELLEN TO 3*)
figure(63);hold on;

% load([pwd,'\R=600P=0.45\GAdat.mat']);
load(fullfile(pwd,'results','v1.mat'));
y=bestofgen{end};
co =[184 77 157]./255;
textz = 'mode 3';
%uncomment if not using matrix file containing correct tunnel tip=3
res =CA_Functions2(y,length(y),432/2,2,0,1,600,.3,3);  %probs,numants,numits*10000,width,infEnergy

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

%% PAPER QUALITY GINI VS. GENERATIONS
figure(55); hold on;
fz= 20;
cc=get(gca,'colororder');
xlab=xlabel('Generation');ylab=ylabel('Gini Coefficient');
gens=0:50;

%equal seed
load([pwd,'\PAPERFIG_GA_R=600_P=0.45 data\GAdatEQ.mat']);
% gens=0:length(bestofgenOUT)-1;

eqGini=zeros(1,length(gens));
for i = 1:length(gens)
    [eqGini(i),~]=Gini(bestofgenOUT{i});
end


%unequal seed
load([pwd,'\PAPERFIG_GA_R=600_P=0.45 data\GAdatUNEQ.mat']);
uneqGini=zeros(1,length(gens));
for i = 1:length(gens)
    [uneqGini(i),~]=Gini(bestofgenOUT{i});
end


%random seed
% load([pwd,'\PAPERFIG_GA_R=600_P=0.45 data\GAdatRAND.mat']);
load(fullfile(pwd,'results','v1.mat'));
randGini=zeros(1,length(gens));
for i = 1:length(gens)
    [randGini(i),~]=Gini(bestofgenOUT{i});
end



plot(gens,eqGini     ,'-','linewidth',3,'color',cc(1,:));
plot(gens,uneqGini   ,'-','linewidth',3,'color',cc(2,:));
plot(gens,randGini   ,'-','linewidth',3,'color',cc(5,:));
text(9,.2,'Equal'  ,'color',cc(1,:),'fontsize',fz);
text(9,.8,'Unequal','color',cc(2,:),'fontsize',fz);
text(9,.4,'Random' ,'color',cc(5,:),'fontsize',fz);
axis([0 length(gens)-1 0 1]);
% set(gca,
% set(gcf,'resize','off');
set(gca,'box','on','linewidth',2,'xtick',[0, 10,20], 'xticklabel',{'0','','20'}, 'ytick',[ .5  1],'yticklabel', {'','1'},'fontsize',fz);
set(xlab,'units','normalized','Position',[.5 -.02 0]);
set(ylab,'units','normalized','Position',[0 .5 0]);
% set(xlab,'position',get(xlab,'position')+[25 .2  0],'margin',1);
axis square
% set(xlab,'margin',3,'clipping','on','position',get(xlab,'position')+[25 .01  0]);

%% GINI cum fraction]
fz=20;
figure(123)
set(gca,'box','on','linewidth',2,'xtick',[0 0.5 1], 'xticklabel',{'0','','1'}, 'ytick',[0 0.5 1],'yticklabel', {'','','1'},'fontsize',fz);
axis square
figure(69); hold on;
load(fullfile(pwd,'PAPERFIG_GA_R=600_P=0.45 data','GAdatUNEQ.mat'));
% y=bestofgen{end};
% %uncomment if not using matrix file containing correct tunnel tip=3
% res =CA_Functions2(y,length(y),432/2,2,0,1,600,.45,3);  %probs,numants,numits*10000,width,infEnergy

z = zeros(1,29);
z(30)=1;
load(fullfile(pwd,'PAPERFIG_GA_R=600_P=0.45 data','GAdatEQ.mat'));
y=bestofgen{1};
res =CA_Functions2(y,length(y),432/2,2,0,1,600,.45,3);  %probs,numants,numits*10000,width,infEnergy
[~,eqLine]=Gini(sum(res.markMatr(:,2:end)));

load(fullfile(pwd,'PAPERFIG_GA_R=600_P=0.45 data','GAdatUNEQ.mat'));
y=bestofgen{1};
res =CA_Functions2(y,length(y),432/2,2,0,1,600,.45,3);  %probs,numants,numits*10000,width,infEnergy
[~,uneqLine]=Gini(sum(res.markMatr(:,2:end)));

% load(fullfile(pwd,'PAPERFIG_GA_R=600_P=0.45 data','GAdatRAND.mat'));
% y=bestofgen{end};
% res =CA_Functions2(y,length(y),432/2,2,0,1,600,.34,3);  %probs,numants,numits*10000,width,infEnergy
% [~,theoryLine]=Gini(sum(res.markMatr(:,2:end)));
cc=get(gca,'colororder');
ucol=cc(1,:);
ecol=cc(2,:);
tcol=cc(5,:);

plot(eqLine(:,1),eqLine(:,2),'linewidth',3,'color',ecol);
plot(uneqLine(:,1),uneqLine(:,2),'linewidth',3,'color',ucol);
% plot(theoryLine(:,1),theoryLine(:,2),'linewidth',3,'color',tcol);

text(.1,.3,'Equal','fontsize',fz,'color',ecol);
text(.4,.2,'Theory','fontsize',fz,'color',tcol);
text(.7,.1,'Unequal','fontsize',fz,'color',ucol);

xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}');
set(gca,'box','on','linewidth',2,'xtick',[0 .5 1], 'ytick',[.5 1],'fontsize',fz);
axis([0 1 0 1]);



%%
load('D:\Dropbox\Social Digging Work\ForPaper\AntExperimentsRepository\Data\AllMarkovMat.mat');
inds=12/240;
%3hrs
inds=round(3/inds);

c=struct;
for i=1:length(ResLorenzN0)
c(i).col=ResLorenzN0(i).col;
c(i).moist=ResLorenzN0(i).moist;
c(i).ants=nnz(sum(ResLorenzN0(i).markmatr1(1:inds,:),1));
end

%% Tunnel length vs. time (*PAPER QUALITY SET TUNNELLEN TO 3*) WITH TIME SCALING TOO
figure(4); hold on;
perc=10;
textS=35;
hr3=3600/.5*3;
min15=60*15/.5;
%exp
%%%%%%%%%%%%%%%mode 1%%%%%%%%%%%%%%%%%%%%%%
tuntip =10;
% load([pwd,'\R=600P=0.45\GAdat.mat']);
load([pwd,'\results\v1.mat']);
y=ones(1,30);
% textz = 'mode 1';
% co =[1 0 0];
comode1 =[1 0 0];
% 432/2 ---> 221
% put pell2grow 200->100!
%
res =CA_Functions2(y,length(y),432/2,2,0,1,600,.4,tuntip);  %probs,numants,numits*10000,width,infEnergy
% [ggruns(i),xy{i}]=Gini(sum(res.markMatr(:,2:end)));
xt=res.tunLength(:,1)*.5/3600;
pellet2grow=res.pell2grow;
pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
pellsTot = sum(res.markMatr(:,2:end),2)/res.pause2dig/pellPerCm;
yt=cumsum(pellsTot);
dat=res.markMatr(1:hr3,2:31);
dat=reshape(dat,min15,30,12);
dat=sum(dat,1);
for(i=1:size(dat,3))
   navgM1(i)=nnz(dat(:,:,i)); 
end
navgM1=mean(navgM1);

% xt2=decimate(xt,100);yt2=decimate(yt-yt(1),100);
mode1x=decimate(xt,100);mode1y=decimate(yt-yt(1),100);
% plot(xt2,yt2,'-','linewidth',4,'color',co,'markersize',15);
% text(2.2,.5,textz,'fontsize',textS,'color',co);

%%%%%%%%%%% mode 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% mode 3 %%%%%%%%%%%%%%%%%%%%%%
%unequal distrib
y=bestofgen{end};
% clear s;
% s.x=1;
% s= SetProbabilities(pwd,30,s,1,0); %9 - probabilities!!! from averaged Lorenz curve for 12 hours
% y=s.prob;
% y=y';
% co =[184 77 157]./255;
% textz = 'mode 3';
comode3 =[184 77 157]./255;
textzmode3 = 'mode 3';

% 432/2 ---> 220
res =CA_Functions2(y,length(y),221,2,0,1,600,.45,tuntip);  %probs,numants,numits*10000,width,infEnergy
% [ggruns(i),xy{i}]=Gini(sum(res.markMatr(:,2:end)));
pellet2grow=res.pell2grow;

pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
pellsTot = sum(res.markMatr(:,2:end),2)/res.pause2dig/pellPerCm;
yt=cumsum(pellsTot);
xt=res.tunLength(:,1)*.5/3600;

dat=res.markMatr(1:hr3,2:31);
dat=reshape(dat,min15,30,12);
dat=sum(dat,1);
for(i=1:size(dat,3))
   navgM3(i)=nnz(dat(:,:,i)); 
end
navgM3=mean(navgM3);

% navgM3=nnz(sum(res.markMatr(1:hr3,2:31),1))

mode3x=decimate(xt,100);mode3y=decimate(yt-yt(1),100)*navgM1/navgM3;

% plot(xt2,yt2,'-','linewidth',4,'color',co,'markersize',15);
% text(1,1.8,textz,'fontsize',textS,'color',co);


%%%%%%%%%%%%%%%%%%mode 3%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%ant expt %%%%%%%%%%%%%%%%%%%%%

%each row represent 3 mins, so for 3 hours of data it will be row 60
if(perc ==10)
    load(['CN16_210um_10per_LONG_NoDupl_Statistics.mat']);
else
    load(['CN16_210um_1per_LONG_EXP2_NoDupl_Statistics.mat']);
end
dat=markmatrYN(1:60,2:31);
dat(dat==-1)=0;
dat=reshape(dat,5,30,12)
dat=sum(dat,1);
for(i=1:size(dat,3))
   c16(i)=nnz(dat(:,:,i)); 
end

if(perc ==10)
    load(['CN12_210um_10per_Converted_NoDupl_TROUBLE_Statistics.mat']);
else
    load(['CN12_210um_1per_LONG_NoDupl_Statistics.mat']);
end

dat=markmatrYN(1:60,2:31);
dat(dat==-1)=0;
dat=reshape(dat,5,30,12)
dat=sum(dat,1);
for(i=1:size(dat,3))
   c12(i)=nnz(dat(:,:,i)); 
end

if(perc ==10)
    load(['CN3_210um_10per_EXP2_NoDupl_Statistics.mat']);
else
    load(['CN3_210um_1per_LONG_NoDupl_Statistics.mat']);
end

dat=markmatrYN(1:60,2:31);
dat(dat==-1)=0;
dat=reshape(dat,5,30,12)
dat=sum(dat,1);
for(i=1:size(dat,3))
   c3(i)=nnz(dat(:,:,i)); 
end



load([pwd,'\LvsT.mat']);
textz='exp';
figure(4); hold on;
% xlabel('{\itT} (h)'); ylabel('{\itL} (cm)');
%lvst = colony/moisture/time/h    time in ? and h in cm
leg = {};
colonies = unique(LvsT(:,1:2),'rows');
% colonies = [16 10];

colonies = colonies(colonies(:,2)==perc,:);
coldig1=[nnz(c3),nnz(c12),nnz(c16)];
% coldig1=[24,23,10];
% coldig10=[14,16,14];
coldig10=[nnz(c3),nnz(c12),nnz(c16)];

if(perc==10)
    colcol=coldig10;
else
    colcol=coldig1;
end
rr=LvsT(LvsT(:,2)==perc,:);

% LVST[colony moist time length]
for i=1:size(colonies,1)
    colony = colonies(i,1);
    
%       rows1=LvsT(:,colony(1)==LvsT(:,1));
%         [rows ~]=find(colony==[]); %all rows of certain colony
%     rows=strmatch(colonies(i,:),LvsT(:,1:2));

% rows=strmatch(colonies(i,:),LvsT(:,1:2));  
     rows=rr(rr(:,1)==colony,:);
     smallTime=rows(rows(:,3)<3.09,:);
     ydig=(smallTime(:,4)-smallTime(1,4))*navgM1/colcol(i);
%      plot(smallTime(:,3),ydig,'.','markersize',35,'linewidth',20);
     
     digs(:,i)=ydig;
% % % % %     %     leg=cat(1,leg,['Col ', num2str(colonies(i,1)),' ',num2str(colonies(i,2)),'%']);
% % % % %     smallTime=find(LvsT(rows,4)<3);
% % % % %     figure;
% % % % %     plot(LvsT(smallTime,3),LvsT(smallTime,4)-LvsT((1),4),'k.','markersize',35,'linewidth',20);

size(smallTime,1)
end
errorbar(smallTime(:,3),mean(digs,2),std(digs,1,2),'-k.','linewidth',2,'markersize',35);
text(2,1.3,textz,'fontsize',textS,'color','k');
% legend(leg);



% xlabel('time (h)'); ylabel('Tunnel Length (cm)');
% title(['Tunnel Length vs. Time N=', num2str(res.numants)]);
% xt=res.tunLength(:,1)*.5/3600;
% yt=res.tunLength(:,2)/2; %bw = .5 cm




plot(mode3x,mode3y,'-','linewidth',4,'color',comode3,'markersize',15);
plot(mode1x,mode1y,'-','linewidth',4,'color',comode1,'markersize',15);

set(gca,'box','on','linewidth',2,'xtick',[0 1 2 3],'ytick',[1 2 3 4],'yticklabel',[],'xticklabel',[],'fontsize',30);
%x=1.793 y=1
set(gcf,'Position',[100,100, 930,600])
% set(gca,'height',.83);
axis([0 3 0 5]);
%%%%%%%%%%%%%%%%%%ant expt %%%%%%%%%%%%%%%%%%%%%
 %% Tunnel length vs. time (*PAPER QUALITY SET TUNNELLEN TO 3*)
% figure(4); hold on;
% 
% textS=35;
% 
% %exp
% %%%%%%%%%%%%%%%mode 1%%%%%%%%%%%%%%%%%%%%%%
% tuntip =3;
% load([pwd,'\R=600P=0.45\GAdat.mat']);
% y=ones(1,30);
% % textz = 'mode 1';
% % co =[1 0 0];
% comode1 =[1 0 0];
% % 432/2 ---> 221
% %
% res =CA_Functions2(y,length(y),432/2,2,0,1,600,.34,tuntip);  %probs,numants,numits*10000,width,infEnergy
% % [ggruns(i),xy{i}]=Gini(sum(res.markMatr(:,2:end)));
% xt=res.tunLength(:,1)*.5/3600;
% pellet2grow=res.pell2grow;
% pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
% pellsTot = sum(res.markMatr(:,2:end),2)/res.pause2dig/pellPerCm;
% yt=cumsum(pellsTot);
% navgM1=nnz(sum(res.markMatr,1));
% % xt2=decimate(xt,100);yt2=decimate(yt-yt(1),100);
% mode1x=decimate(xt,100);mode1y=decimate(yt-yt(1),100);
% % plot(xt2,yt2,'-','linewidth',4,'color',co,'markersize',15);
% % text(2.2,.5,textz,'fontsize',textS,'color',co);
% 
% %%%%%%%%%%% mode 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%% mode 3 %%%%%%%%%%%%%%%%%%%%%%
% %unequal distrib
% y=bestofgen{end};
% % clear s;
% % s.x=1;
% % s= SetProbabilities(pwd,30,s,1,0); %9 - probabilities!!! from averaged Lorenz curve for 12 hours
% % y=s.prob;
% % y=y';
% % co =[184 77 157]./255;
% % textz = 'mode 3';
% comode3 =[184 77 157]./255;
% textzmode3 = 'mode 3';
% 
% % 432/2 ---> 220
% res =CA_Functions2(y,length(y),221,2,0,1,600,.45,tuntip);  %probs,numants,numits*10000,width,infEnergy
% % [ggruns(i),xy{i}]=Gini(sum(res.markMatr(:,2:end)));
% pellet2grow=res.pell2grow;
% 
% pellPerCm = pellet2grow*2*2; %tunnwidth=2,another 2 for 1site=.5 cm
% pellsTot = sum(res.markMatr(:,2:end),2)/res.pause2dig/pellPerCm;
% yt=cumsum(pellsTot);
% xt=res.tunLength(:,1)*.5/3600;
% NANTSM3=nnz(sum(res.markMatr,1));
% mode3x=decimate(xt,100);mode3y=decimate(yt-yt(1),100)*navgM1/navgM3;
% 
% % plot(xt2,yt2,'-','linewidth',4,'color',co,'markersize',15);
% % text(1,1.8,textz,'fontsize',textS,'color',co);
% 
% %%%%%%%%%%%%%%%%%%mode 3%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% load([pwd,'\LvsT.mat']);
% textz='exp';
% figure(4); hold on;
% % xlabel('{\itT} (h)'); ylabel('{\itL} (cm)');
% %lvst = colony/moisture/time/h    time in ? and h in cm
% leg = {};
% colonies = unique(LvsT(:,1:2),'rows');
% % colonies = [16 10];
% perc= 10;
% colonies = colonies(colonies(:,2)==perc,:);
% coldig1=[24,23,10];
% coldig10=[14,16,14];
% 
% rr=LvsT(LvsT(:,2)==perc,:);
% 
% % LVST[colony moist time length]
% for i=1:size(colonies,1)
%     colony = colonies(i,1);
%     
% %       rows1=LvsT(:,colony(1)==LvsT(:,1));
% %         [rows ~]=find(colony==[]); %all rows of certain colony
% %     rows=strmatch(colonies(i,:),LvsT(:,1:2));
% 
% % rows=strmatch(colonies(i,:),LvsT(:,1:2));  
%      rows=rr(rr(:,1)==colony,:);
%      smallTime=rows(rows(:,3)<3.09,:);
%      ydig=(smallTime(:,4)-smallTime(1,4))*navgM1/coldig10(i);
% %      plot(smallTime(:,3),ydig,'.','markersize',35,'linewidth',20);
%      
%      digs(:,i)=ydig;
% % % % % %     %     leg=cat(1,leg,['Col ', num2str(colonies(i,1)),' ',num2str(colonies(i,2)),'%']);
% % % % % %     smallTime=find(LvsT(rows,4)<3);
% % % % % %     figure;
% % % % % %     plot(LvsT(smallTime,3),LvsT(smallTime,4)-LvsT((1),4),'k.','markersize',35,'linewidth',20);
% 
% size(smallTime,1)
% end
% errorbar(smallTime(:,3),mean(digs,2),std(digs,1,2),'-k.','linewidth',2,'markersize',35);
% text(2,1.3,textz,'fontsize',textS,'color','k');
% % legend(leg);
% 
% 
% 
% % xlabel('time (h)'); ylabel('Tunnel Length (cm)');
% % title(['Tunnel Length vs. Time N=', num2str(res.numants)]);
% % xt=res.tunLength(:,1)*.5/3600;
% % yt=res.tunLength(:,2)/2; %bw = .5 cm
% 
% 
% 
% 
% plot(mode3x,mode3y,'-','linewidth',4,'color',comode3,'markersize',15);
% plot(mode1x,mode1y,'-','linewidth',4,'color',comode1,'markersize',15);
% 
% set(gca,'box','on','linewidth',2,'xtick',[0 1 2 3],'ytick',[1 2 3 4],'yticklabel',[],'xticklabel',[],'fontsize',30);
% %x=1.793 y=1
% set(gcf,'Position',[100,100, 930,600])
% % set(gca,'height',.83);
% axis([0 3 0 5]);

%% plot growth rate vs ants with both modes (*PAPER QUALITY*) tunneltip=10?
clear all;
outer=jet(9);
figure(11); hold on;
% titl= 'Tunnel Growth Rate vs. N';
xlabel('{\itn}'); ylabel('{\itV} (cm/h)');
p='0.3'
for eq=[0 1] %% zero for inequal work
    xText=60;
    if(eq)
        
        fold = ['\p',p,'eqpell200\'];
        %     fold ='\eqResSavesInfEnergyPell600\';
        %     fold ='\eqResSavesRegEnergyPell200\';
        textz = 'Equal';
        yText = .25;
        mark = 's';
        co = [1 0 0];
    else
        fold =['\p',p,'uneqpell200\'];
        %     fold ='\uneqResSavesInfEnergyPell600\';
        %     fold ='\uneqResSavesRegEnergyPell200\';
        mark = 'o';
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
    plot(antNum,V,mark,'LineWidth',2,'color',co,'markersize',8,'linewidth',4);
    
    %     plot(antNum,V,[mark, '-'],'LineWidth',2,'color',outer(str2double(p)*10,:),'MarkerFaceColor',co,'markersize',6,'linewidth',1);
    
    axis([0 100 0 0.7]);
    text(xText,yText,textz,'fontsize',35,'color',co); %%%UNCOMMENT%%%%
    set(gca,'box','on','linewidth',2,'xtick',[0 50 100],'ytick',[.3 .6],'fontsize',35);
    xlabh=get(gca,'xlabel');
    % set(xlabh,'Position',get(xlabh,'Position')+[0 .01 0]);
end

% x=2.6 y = 1
set(gcf,'Position',[100,100, 1300/1.2, 500/1.2])

%% plot q_bar vs rho_bar ****PAPER QUALITY*************SET tunneltip=10
figure(14); hold on;
titl='Traffic (CA)';
p='0.35'
% h1=xlabel('$\bar{\rho}$ (ants/site)');
% h2=ylabel('$\bar{q}$ (ants/min)');
fsize=24;
% h1=xlabel('$\bar{\rho}$ (ants/BL)');
% h2=ylabel('$\bar{q}$ ($\frac{ants}{BL{\cdot}min}$)');
% h2=ylabel('$\bar{q}$ (ants/(BL$\cdot$min))','fontsize',fsize);

% set(h1,'interpreter','Latex','FontWeight','bold');
% set(h2,'interpreter','Latex','FontWeight','bold','fontsize',fsize);
% eq=0; %% zero for inequal work
for eq = [0 1]
    
    if(eq)
        %     fold ='\eqResSavesInfEnergy\';
        %     fold ='\eqResSavesInfEnergyPell600\';
        fold = ['\p',p,'eqpell200\'];
        mark = 's';
        lw=2;
    else
        %     fold ='\uneqResSavesInfEnergy\';
        %     fold ='\uneqResSavesInfEnergyPell600\';
        fold =['\p',p,'uneqpell200\'];
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
    cmap = jet(numcolors);
    stopTime =2.75;
    startTime=.25;
    minz = 60*(stopTime-startTime); %%%uncomment
    %   minz = 60*(stopTime-startTime);
    stopTime = stopTime*7200;
    startTime=startTime*7200;
    outer=jet(9);
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
        
        
        a=sum(floor(sum(res.markMatr(startTime:stopTime,2:end))/res.pause2dig));
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
        plot(rhom,qm,mark,'MarkerSize',15,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k','LineWidth',lw);
        %         plot(rhom,qm,mark,'MarkerSize',6,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',outer(str2double(p)*10,:),'LineWidth',1);
        
    end
    haxis = gca;
    set(gca,'fontsize',55);
    figText(gcf,16)
    colormap(cmap);
    
%     caxis([1 numcolors])
    
%     cbarHandle = colorbar('YTick',...
%         1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors,...
%         'YTickLabel',ants, 'YLim', [1 numcolors]); %[]=ants
    
    % set(gca,'XTick',[0.25 0.5],'YTick',[0 .1 .2],'XLim',[0.045 0.55],'YLim',[0,0.2]);
    set(gca,'XTick',[0.15 0.3],'YTick',[0.025 .05],'yTicklabel',[], 'xTicklabel',[],'fontsize',fsize);
%     set(cbarHandle,'fontsize',18');
    axis([0 .4 0 .075])
    if(res.infEnergy)
        En=' E=\infty';
    else
        En=[];
    end
    % title([titl,' (',num2str(stopTime/7200),'h)',En]);
    % title([titl,En]);
    set(gca,'box','on','linewidth',2);
    set(gcf,'position',[100,100,426,429]);
end

%% ANT DATA WITH THEORY LINE
figure(6);hold on;


fz = 24;
clear xy
load([pwd,'\R=600P=0.45\GAdat.mat']);

runs = 5;
y=bestofgen{end};
for(i=1:runs)
    res =CA_Functions2(y,length(y),432/2,2,0,1,600,.45,10);  %probs,numants,numits*10000,width,infEnergy
    [ggruns(i),xy{i}]=Gini(sum(res.markMatr(:,2:end)));
    %      CA_FunctionsWill(prob,length(prob),ni,2,0,1,rec,ptt);
end
figure(6);
yy=cell2mat(cellfun(@(x) x(:,2),xy,'UniformOutput',0));
meanY=mean(yy,2);
% meanY=mean([a{:}],2);
meanG=mean(ggruns);
err=std(yy,1,2);
cc=parula(15);
shadedErrorBar(xy{1}(:,1),meanY,err,{'Color',cc(10,:),'LineWidth',3},1);
% text(.651,.829,['Theory'],'Color',cc(10,:),'fontsize',fz);
% xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}');
% figText(gcf,18);
% xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}');
set(gca,'box','on','linewidth',2);
set(gcf,'Position',[1 1 594 562]);



for i=1:2 %1 =1% 2=10%
    if i==1
        load('AntData\antData1p.mat');
    else
        load('AntData\antData10p.mat');
    end
    if per==1
        %    text(.65,.1,['{\itW} = 0.01 G=',num2str(mwG)],'Color','red','FontSize',18)
        %         text(.419,.555,['{\itW} = 0.01'],'Color','red','FontSize',fz)
        ccl = [1,0,0];
    else
        cc=parula(13);
        ccl=cc(2,:);
        cct=cc(1,:);
        %     text(.5,.5,['{\itW} = 0.1 G=',num2str(mwG)],'Color',cct,'FontSize',18)
        %         text(.7,.1,['{\itW} = 0.1'],'Color',cct,'FontSize',fz)
        lcol = 'b';
    end
    for i=1:expNum
        if strcmp(h20,'1 and 10%')
            
        else
            [gg(i),xy{i}]=Gini(ResLorenzN0(i).totalMarks(end,:));
            
            xy{i}(end,:)=[1,1];
            xyInterp(:,i)=interp1(xy{i}(:,1),xy{i}(:,2),linspace(0,1,30));
        end
    end
    mw = mean(xyInterp,2);
    [mwG,~]=Gini(diff(mw));
    meanG=mean(gg);
    % xxy=(xy(:,2,:)');
    errw= [std(xyInterp,0,2)];
    shadedErrorBar(linspace(0,1,30),mw,errw,{'-','Color',ccl,'LineWidth',1.2},1);
    set(gca,'YTick',[0.5 1],'XTick',[0 0.5 1],'yTicklabel',[],'XTicklabel',[],'XLim',[0 1],'YLim',[0 1],'fontsize',fz);
    
end
% plot theory line

% w=1027.0 h=595;
%% figure with inset optimal lorenz with gini vs gen inset
figure(69); hold on;
load([pwd,'\R=600P=0.45\GAdat.mat']);
y=bestofgen{end};
%uncomment if not using matrix file containing correct tunnel tip=3
res =CA_Functions2(y,length(y),432/2,2,0,1,600,.45,3);  %probs,numants,numits*10000,width,infEnergy

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
text(.1,.3,'Equal','fontsize',fz,'color',ecol);
text(.4,.2,'Theory','fontsize',fz,'color',tcol);
text(.7,.1,'Unequal','fontsize',fz,'color',ucol);

xlabel('\it{Cum. fraction of workers}');ylabel('{\itCum. fraction of work}');
set(gca,'box','on','linewidth',2,'xtick',[0 .5 1], 'ytick',[.5 1],'fontsize',fz);
axis([0 1 0 1]);

% INSET PART GINI VS GEN BELOW

axes('position',[.2 .6 .3 .3]);
hold on;
cc=get(gca,'colororder');
xlab=xlabel('Generation');ylab=ylabel('Gini');

%equal seed
load([pwd,'\PAPERFIG_GA_R=600_P=0.45 data\GAdatEQ.mat']);
gens=0:20;
eqGini=zeros(1,length(gens));
for i = 1:length(gens)
    [eqGini(i),~]=Gini(bestofgenOUT{i});
end


%unequal seed
load([pwd,'\PAPERFIG_GA_R=600_P=0.45 data\GAdatUNEQ.mat']);
uneqGini=zeros(1,length(gens));
for i = 1:length(gens)
    [uneqGini(i),~]=Gini(bestofgenOUT{i});
end


%random seed
load([pwd,'\PAPERFIG_GA_R=600_P=0.45 data\GAdatRAND.mat']);
randGini=zeros(1,length(gens));
for i = 1:length(gens)
    [randGini(i),~]=Gini(bestofgenOUT{i});
end



plot(gens,eqGini     ,'-','linewidth',2,'color',ecol);
plot(gens,uneqGini   ,'-','linewidth',2,'color',ucol);
plot(gens,randGini   ,'-','linewidth',2,'color',tcol);
text(9,.15,'Equal'  ,'color',ecol,'fontsize',15);
text(9,.75,'Unequal','color',ucol,'fontsize',15);
text(9,.4,'Random' ,'color',tcol,'fontsize',15);
axis([0 length(gens)-1 0 1]);
% set(gca,
set(gcf,'resize','off');

set(xlab,'units','normalized','Position',[.5 -.02 0]);
set(ylab,'units','normalized','Position',[0 .5 0]);
set(gca,'box','on','linewidth',2,'xtick',[0, length(gens)-1], 'ytick',[1],'fontsize',15);
set(gcf,'Position',[1,1, 800, 600])


%% robot 1/n*DE/DN vs. n,
%[mode1 mode2 mode3]=[r b m]
%mode 1=[active, (very small reversal necessary to not have permajam)]
%mode 2=[active, reversal]
%mode 3=[lazy, no reversal]
%mode 4=[lazy, reversal]


%mode 1:

% openfig('A:\GA\dEdNn.fig','reuse')
figure (1);
openfig('A:\GA\N.fig','reuse');
% openfig('A:\GA\dNdT.fig','reuse')
%[probs, numants, numits*100, tunwidth, infEnergy, placeHolder(1), rechsteps, prob2turn, tuntip]
for modez=3
    co =[1 0 0;0,0,1;184/255 77/255 157/255];
    dtt=30;
    st = 9; %simtime
    % probs =[50^dtt,50^dtt,50^dtt,80^dtt];
    % probs =[0.285,0.014,0.888,0.569];
    probs =[0.883,0.006,0.029,0.930];
    robs=[1:4];
    it = 1;
    % kk =1;
    jIts=1;
    pturn =[0.001, 0.01, 0.001];
    for j = 1:4
        tic
        pts('ants:',j);
        for it = 1:50
            if modez>1 %if lazy mode
                %                 y=probs(1:j)/sum(probs(1:j));
                y=probs(1:j);
            else
                y=ones(1,j);
            end
            [res(it) ~]=RobotCA(modez,j);
            de(it)=abs(sum(res(it).indEnergy(1,:))-sum(res(it).indEnergy(end,:))); %Efinal-Einitial
            dN(it)=sum(res(it).pellets(end,:));
            dedNn(it)=(1/j)*(de(it)/dN(it));
            %         find(diff(res(it).indEnergy)>0)
        end
        removeAmt = 5;
        [~,indz]=sort(dedNn,'descend');
        nindz=[indz(1:removeAmt) indz(end-removeAmt:end)];
        dedNn(nindz)=[];%remove largest and lowest 5
        dedNn_m(jIts) = mean(dedNn);
        dedNn_err(jIts)=std(dedNn);
        
        [~,indz]=sort(dN,'descend');
        nindz=[indz(1:removeAmt) indz(end-removeAmt:end)];
        dN(nindz)=[];
        N_m(jIts) = mean(dN);
        N_err(jIts) = std(dN);
        
        jIts=jIts+1;
        toc;
    end
    % figure(1);
    % hold on;
    %
    % % set(gcf,'resize','off');
    %
    % xlabel('Number of Robots');
    % ylabel('$\displaystyle\frac{1}{n}\frac{dE}{dN}$','interpreter','latex');
    %
    % % errorbar(robs,dedNn_m,dedNn_err,'color',co(modez,:),'linewidth',2);
    % errorbar(robs,dedNn_m,dedNn_err,'linewidth',2,'color',co(modez,:));
    % figText(gcf,14,1);
    
    
    figure(2);
    hold on;
    
    xlabel('Number of Robots');
    ylabel('$\displaystyle N, deposits$','interpreter','latex');
    errorbar(robs,N_m,N_err,'linewidth',2,'color',co(modez,:));
    
    % figure(3)
    % hold on;
    % xlabel('Number of Robots');
    % ylabel('$\displaystyle\frac{dN}{dt}, \frac{deposits}{minute}$','interpreter','latex');
    % dndt= 30/900;
    % errorbar(robs,dndt*N_m,N_err*dndt,'linewidth',2,'color',co(modez,:));
    
    clear all
end