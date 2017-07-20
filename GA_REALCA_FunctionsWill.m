% function [ResEnergy]=CA_Functions(prob,numantsALL)
function [excavated]=GA_REALCA_FunctionsWill(proba,tt,rech,ptt,tuntip)%%NORMAL
% function [excavated]=GA_REALCA_FunctionsWill(proba,tt %%RECH AND PTT
% tic
rng('shuffle');
DRAW = 0;
movieDir='A:\movie';
% numantsALL=[2,3,5,8,10,12,15,18,20,25,30];

tw=2;
energyMult=1;

%%%%%RECH AND PTT%%%%%%%%
% recharge_steps=proba(31);
% prob_turn=proba(32);
% numants=length(proba)-2;

%%%NORMAL%%%%%%%
recharge_steps=rech; %orginally =20steps now 20 mins
prob_turn=ptt; %has to be 0.34
numants = length(proba);


cmap=[1 0 0; 0.8 0.8 0; 1 1 1; 1 0.5 0;rand(numants,3)];

% clearvars -except numants DRAW run workfolder schRes ResEnergy groupEnergy...
%     indEnergy tunLength numantsALL z proba tt tw 1 energyMult cmap rech
steps2drop=20; % # of steps to drop the pellet
pause2dig=25; %number of steps ant spends digging - 25



% important parameter - defines how much time the colony actually is actually digging


pellet2grow=200; %100/tunnelsize;% should be 100/tunnel size - number of pellets to increase the tunnel tip


exc_ECost=3;
Transp_ECost=4;%SHOULD BE 4!!!
walk_ECost=1;

tunnelsize = tw;
tunneltip=tuntip;
roadlength = 100;
countSpot = roadlength-2;
iterations = tt*100; % .5 secs per iteration 345600 = 48 hours

elow_limit=(walk_ECost*roadlength+Transp_ECost*roadlength+exc_ECost*pause2dig)/2;%*2
etop_limit=elow_limit*10 ;

recharge=round((etop_limit-elow_limit)/recharge_steps);% energy gained after 1 recharging step
pelletcount=0; % number of pellets excavated in total

coeffP=1;
%         coeffP=0.93; %=within first 12 hours ants dig in the tunnel only 93% (48h-68%) of the time -93 for future: ants are digging coeffP*prob2dig
%0 - unexcavated, 1 - excavated, 2 and -2 - ants
penalty=0;
% all1=0; % if all 1 is set to 1, ants are digging continuously, if to 0 - with probabilities P
uselessRuns=0;% runs which did not result in the pellet excavation
% set initial conditions
prob_lateral = 0.52;
prob_forward = 1;


%prob to move with pellet
if penalty==1
    prob_lateralp = 0.7;
    prob_forwardp = 0.4;
else
    prob_lateralp = 0.52; % was 0.4 with penalty NO DIFFERENCE FOR NOW
    prob_forwardp = 1; % was 0.7
end

%%%%%%%initialize road%%%%%%%%%%%%
%       road = Initialize(tunnelsize, tunneltip, roadlength, numants);
road = zeros(tunnelsize, roadlength); %0-sand
road(1:tunnelsize,roadlength-tunneltip:roadlength)=1; %1 - excavated road
%%%%%%%%%%%%%%%%%%%%%

%%%this was set on at 7/20/17
%         if(DRAW)
%             
%             imagesc(road);
%             drawnow;
%         end

%initially ALL ants have to be set with coordinates equal to the
%tunnel length

systemState.x(1:numants,1)=roadlength;  %x of ants
systemState.y(1:numants,1)=randi([1,tunnelsize],1,numants); % of ants
systemState.digpause=zeros(numants,1);  %contains all ants, which present at 3 BL away from the tunnel face          %3 - pause for digging
systemState.pellet=zeros(numants,1);                 %4 - 1 for pellet/0 for no pellet
systemState.energy(1:numants,1)=etop_limit*energyMult;          %5 - energy
systemState.globalindx(1:numants,1)=1:numants;       %6 - indexes in global system
systemState.restflag=ones(numants,1);                %7 - Flag to stop digging (1/ if ant RESTs/recharges energy, 0 when ant walking): FOR START EVERYONE IS NOT AT THE TUNNEL
systemState.time2drop=zeros(numants,1);              %8 - # of steps ant spent droping the pellet
systemState.prob=zeros(numants,1);
systemState.atFace=zeros(numants,1);
resting=zeros(numants,1);
%         systemState.atFace=zeros(numants,1);
%         systemState.prob = prob;


%         systemState.prob=p.^(1.75);
systemState.prob=proba(1:30);
systemState.prob=systemState.prob./sum(systemState.prob);

% counter for total energy
m=1;
%         growth_energy = zeros(m,3);
% growth_energy(m,1:3)=0;
growth_energy=zeros(iterations*numants,3);
%%%%%%%%%%%%%%
%added code which preallocates variables which change size on each
%iteration previously
markmatr = zeros(iterations,numants);
markmatrN0 = zeros(iterations,numants);
tunnel_length = zeros(iterations,2);
tunnel_length(:,1)= 1:iterations;
energy = zeros(iterations,numants);
pellet = zeros(iterations,numants);
probs = zeros(iterations,numants);
atFace = zeros(iterations,numants);
prevResting = systemState.restflag;
tunTime = zeros(1,numants);
flow = zeros(iterations,2);
occupied = zeros(iterations,2);
%         groupEnergy=zeros(
%%%%%%%%%%%%%%

for kk=1:iterations
    
    ind = randperm(numants);
    systemState.x(:)= systemState.x(ind);
    systemState.y(:)= systemState.y(ind);
    x=systemState.x(:);
    y=systemState.y(:);
    systemState.digpause(:)=systemState.digpause(ind);                %3
    systemState.pellet(:)=systemState.pellet(ind);                    %4
    systemState.energy(:)=systemState.energy(ind);                    %5
    systemState.globalindx(:)=systemState.globalindx(ind);            %6
    systemState.restflag(:)=systemState.restflag(ind);                %7
    systemState.time2drop(:)=systemState.time2drop(ind);              %8
    systemState.prob(:)=systemState.prob(ind);
    systemState.atFace(:)=systemState.atFace(ind);
    prevResting(:)=prevResting(ind);
    resting(:)=resting(ind);
    tunTime(:)=tunTime(ind);
    markmatr(kk,1)=kk;% contains all ants, which present at 3 BL away from the tunnel face
    markmatrN0(kk,1)=kk; % contains all ants, which DIG(!) at 3 BL away from the tunnel face
    moved = zeros(1,2);
    
    %              %%%%%%%%%%%%
    %              vdir= road(sub2ind(size(road),y,x))./abs(road(sub2ind(size(road),y,x)));
    %              rechargeInds = find(systemState.restflag(:,1)==1 &  systemState.pellet(:,1)~=1);
    %              systemState
    %              %%%%%%%%%%%%
    
    for jj=1:numants
        %figure out direction of motion +1 up /-1 down
        vdir = road(y(jj), x(jj))/abs(road(y(jj), x(jj))); %negative
        %%%%%%%%%%%%%%%%%%%recharge%%%%%%%%%%%%%%%%
        if systemState.restflag(jj)
            if systemState.pellet(jj)
                %PELLETDEPOSIT
                if systemState.time2drop(jj,1) < steps2drop
                    systemState.time2drop(jj,1)=systemState.time2drop(jj,1)+1;
                else
                    if systemState.prob(jj)*coeffP>=rand% ant goes to the tunnel
                        if road(y(jj), x(jj))==1
                            road(y(jj), x(jj))=-2;
                            systemState.restflag(jj,1)=0;
                            systemState.pellet(jj,1)=0;
                        else
                            systemState.time2drop(jj,1)=systemState.time2drop(jj,1)-1;
                        end
                    end
                end
                
            else %no pellet
                %RECHARGE
                if systemState.energy(jj,1)>=etop_limit
                    r=rand;
                    if  systemState.prob(jj)*coeffP>=r% ant goes to the tunnel
                        %after ant recharges she starts moving towards tunnel
                        if road(y(jj), x(jj))==1
                            road(y(jj), x(jj))=-2;
                            systemState.restflag(jj,1)=0;
                        end
                    end
                else
                    systemState.energy(jj,1)=systemState.energy(jj,1)+recharge;
                end
            end
        else %%systemState.restflag==0
            if systemState.digpause(jj)~=0
                [systemState,road,tunneltip,growth_energy,m,~,pelletcount,pellet2grow] = Digging(systemState,jj,pause2dig,exc_ECost,walk_ECost,growth_energy,m, road, x, y, roadlength, tunneltip, kk, pelletcount, pellet2grow);
            else %==0
                if vdir<0 && systemState.pellet(jj)==0
                    if x(jj) == (roadlength-tunneltip)
                        %SetToDigTip
                        systemState.energy(jj,1)=systemState.energy(jj,1)+walk_ECost;
                        growth_energy(m,2)=growth_energy(m,2)-walk_ECost;
                        systemState.digpause(jj,1)=1;
                    elseif x(jj)>(roadlength-tunneltip) && ...
                            x(jj)<=(roadlength-tunneltip+2) && ...
                            road(y(jj),x(jj)-1)==-2 && ...
                            ~isempty(find(road(:,x(jj))==1,1))
                        %SetToDigNearTip
                        jjfa=find(systemState.x(:)==x(jj)-1 & systemState.y(:)==y(jj));% the number jj of the front ant to check if she digs
                        
                        if systemState.digpause(jjfa)~=0 && systemState.digpause(jjfa)~=pause2dig %if ant infront digs and not ready to leave the tunnel
                            systemState.energy(jj,1)=systemState.energy(jj,1)+walk_ECost;
                            growth_energy(m,2)=growth_energy(m,2)-walk_ECost;
                            systemState.digpause(jj,1)=1; %stay and dig
                        end
                    end
                end
            end
        end
        
        if systemState.restflag(jj)~=1 && systemState.digpause(jj)==0
            [road, systemState, growth_energy, ~,uselessRuns,moved] = Walking(systemState, jj, prob_lateral, prob_forward, prob_lateralp, prob_forwardp, walk_ECost,growth_energy,m, Transp_ECost, road, x, y, roadlength, tunneltip, elow_limit, prob_turn,tunnelsize,uselessRuns,countSpot,moved);
        end
        %order digging / set to digging/ walking is important
        if(x(jj) == (roadlength-tunneltip))
            systemState.atFace(jj)=systemState.atFace(jj)+1;
        end
        
        if x(jj)>=(roadlength-tunneltip) && x(jj)<=(roadlength-tunneltip+2)
            markmatr(kk,systemState.globalindx(jj)+1)=1;
            if systemState.digpause(jj,1)~=0
                markmatrN0(kk,systemState.globalindx(jj)+1)=1;
            end
        end
        
        if(x(jj)<roadlength)
            tunTime(jj)=tunTime(jj)+1;
        end
    end % end of ant cycle
    
    
    
    if(DRAW && mod(kk,1)==0)
        figure(2);
        newR = road+2;
        for i = 1:numantsALL
            id = find(systemState.globalindx==i);
            if(systemState.digpause(id)>1)
                newR(systemState.y(id),...
                    systemState.x(id))=1;
            else
                newR(systemState.y(id),...
                    systemState.x(id))=i+4;
            end
            
        end
        image(newR);
        colormap(cmap);
        colorbar;
        drawnow;
        saveas(gcf,fullfile(movieDir,[num2str(kk),'.png']));
        
    end
    
    %             tunnel_length(kk,1)=kk;
    tunnel_length(kk,2)=tunneltip;
    %change in energy of every single ant on each cycle step
    energy(kk,systemState.globalindx(:))=systemState.energy(:);
    pellet(kk,systemState.globalindx(:))=systemState.pellet(:);
    probs(kk,systemState.globalindx(:))=systemState.prob(:);
    atFace(kk,systemState.globalindx(:))=systemState.atFace(:);
    flow(kk,:)= moved;
    occupied(kk,:)= [(abs(road(1,countSpot))==2), (abs(road(2,countSpot))==2)];
    if tunneltip==1
        break;
    end
    
    %             probabil(systemState.globalindx(:))= systemState.prob(:);
    %             atFace(systemState.globalindx(:))=atFace(:);
    
    
    restInds= find(systemState.restflag-prevResting==-1); %resting=1 to walking=0 0-1=-1
    resting(restInds)=resting(restInds)+1;
    prevResting=systemState.restflag(:);
    
end
excavated=-sum(sum(markmatr(:,2:end)))/1000;
% ResEnergy.groupEnergy=growth_energy(1:size(growth_energy)-1,:); % total energy expenditure per unit of tunnel growth: 1 -tunneltip; 2 -energy expenditure; 3 - number of steps kk
% ResEnergy.numants=numants;
%
% ResEnergy.indEnergy = energy; % energy of individual ants per 1 step of cycle
% ResEnergy.tunLength=tunnel_length; % tunnel growth
% %         markMatr(:,:,run)=markmatr;
% ResEnergy.markMatr=markmatrN0;
% ResEnergy.pellets=[zeros(1,numants); cumsum(diff(pellet)==-1)];
% ResEnergy.pall1=all1;
% ResEnergy.coeffP=coeffP;
% ResEnergy.penalty=penalty;
% ResEnergy.uselessRuns= uselessRuns;
% ResEnergy.numants=numants;
%
% tunTime(systemState.globalindx(:))=tunTime(:);
% ResEnergy.tunTime =tunTime;
%
% ResEnergy.prob=probs;
% ResEnergy.atFace= atFace;
% ResEnergy.resting=resting;
% ResEnergy.flow = flow;
% ResEnergy.occupied = occupied;
% ResEnergy.pell2grow = pellet2grow;
% ResEnergy.infEnergy = energyMult>1;
% ResEnergy.equalDis     = all(probs(1,1)==probs(1,:));


%         ResEnergy(schRes).markMatr=markMatr;   %%? next line it changes?

% toc