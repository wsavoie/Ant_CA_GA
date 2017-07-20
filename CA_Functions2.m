% function [ResEnergy]=CA_Functions(prob,numantsALL)
function [ResEnergy]=CA_Functions2(varargin)
% if nargin<6%so I can supress print out if I want
tic
% end
rng('shuffle');
movieDir='A:\movie';
DRAW = 0;
workfolder=pwd;
% schRes=1;
% numantsALL=[2,3,5,8,10,12,15,18,20,25,30];
numantsALL= 10;
tt=5;
tw=2;
runIts = 1;
energyMult=1;
rech = 600; %default recharge steps
ptt=.7; %default probability to turn
tuntip=20;
if(nargin>1)
    proba = varargin{1};
    numantsALL=varargin{2};
    if(nargin>=3)
        tt=varargin{3};
    end
    if(nargin>=4)
        tw=varargin{4};
    end
    if(nargin>=5)
        energyMult=varargin{5}*10000+1; %1 if {5}=0, 10001 else
    end
    if(nargin>=7)
        rech=varargin{7}; %1 if {5}=0, 10001 else
    end
    if(nargin>=8)
        ptt=varargin{8}; %1 if {5}=0, 10001 else
    end
    if(nargin>=9)
        tuntip=varargin{9}; %1 if {5}=0, 10001 else
    end
end

ResEnergy(size(numantsALL,2))=struct;
for z=1:size(numantsALL,2)
    numants=numantsALL(z);
    cmap=[1 0 0; 0.8 0.8 0; 1 1 1; 1 0.5 0;jet(numantsALL)];
    %     disp(z);
    for run=runIts
        clearvars -except numants DRAW run workfolder schRes ResEnergy groupEnergy...
            indEnergy tunLength numantsALL z proba tt tw runIts energyMult cmap rech ptt tuntip
        countSize= 4;
        steps2drop=20; % # of steps to drop the pellet
        pause2dig=25;%%25; %number of steps ant spends digging - 25
        recharge_steps=rech; %orginally =20steps now 20 mins
        %         recharge_steps=20; % important parameter - defines how much time the colony actually is actually digging
        
        %  pellet2grow=100; %100/tunnelsize;% should be 100/tunnel size - number of pellets to increase the tunnel tip
        %1/30/17
        pellet2grow=200; %100/tunnelsize;% should be 100/tunnel size - number of pellets to increase the tunnel tip
        
        
        exc_ECost=3;
        Transp_ECost=4;%SHOULD BE 4!!!
        walk_ECost=1;
        prob_turn=ptt; %has to be 0.34
        
        tunnelsize = tw;
        tunneltip=tuntip;
        roadlength = 100;
        countSpot = roadlength-2;
        iterations = tt*100; % .5 secs per iteration 345600 = 48 hours
        
        elow_limit=(walk_ECost*roadlength+Transp_ECost*roadlength+exc_ECost*pause2dig)/2;
        etop_limit=elow_limit*10 ;
        
        recharge=round((etop_limit-elow_limit)/recharge_steps);% energy gained after 1 recharging step
        pelletcount=0; % number of pellets excavated in total
        
        coeffP=1;
        %         coeffP=0.93; %=within first 12 hours ants dig in the tunnel only 93% (48h-68%) of the time -93 for future: ants are digging coeffP*prob2dig
        %0 - unexcavated, 1 - excavated, 2 and -2 - ants
        penalty=0;
        all1=0; % if all 1 is set to 1, ants are digging continuously, if to 0 - with probabilities P
        uselessRuns=0;% runs which did not result in the pellet excavation
        % set initial conditions
        prob_lateral =.52;% 0.52;
        prob_forward = 1;
        
        
        %prob to move with pellet
        if penalty==1
            prob_lateralp = 0.7;
            prob_forwardp = 0.4;
        else
            prob_lateralp =prob_lateral;
            %             prob_lateralp = .1;%0.52; % was 0.4 with penalty NO DIFFERENCE FOR NOW
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
        
        if(nargin>1)
            p = proba;
            if all(p(1)==p) %if equal
                systemState.prob = p;
                %                 if nargin<6%so I can supress print out if I want
                pts('ants = ',length(p));
                %                 end
                
            else
                
                systemState.prob=proba;
                systemState.prob=systemState.prob./sum(systemState.prob);
                
                % % %                 p=p/sum(p);
                % %                 systemState.prob=p.^(1.75);
                % %                 systemState.prob=systemState.prob./sum(systemState.prob);
                % % %                 systemState.prob = p;
                %                 if nargin<6
                pts('ants = ',length(p));
                %                 end
                
                %                 pts(systemState.prob);
            end
            
        else
            systemState = SetProbabilities(workfolder,numants,systemState,coeffP,all1); %9 - probabilities!!! from averaged Lorenz curve for 12 hours
        end
        
        % counter for total energy
        m=1;
        %         growth_energy = zeros(m,3);
%         growth_energy(m,1:3)=0; %%initialize to some very large value we won't know final size
        growth_energy=zeros(iterations*numants,3);
        %%%%%%%%%%%%%%
        %added code which preallocates variables which change size on each
        %iteration previously
        markmatr = zeros(iterations,numants+1);
        markmatrN0 = zeros(iterations,numants+1);
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
        density = zeros(iterations,2);
        %         groupEnergy=zeros(
        groupEnergy = zeros(size(growth_energy,1),3,length(runIts));
        indEnergy = zeros(iterations,numants,length(runIts)); % energy of individual ants per 1 step of cycle
        tunLength = zeros(iterations,2,length(runIts)); % tunnel growth
        %         markMatr=zeros(iterations,numants+1,length(runIts));
        %         markMatrN0=zeros(iterations,numants+1,length(runIts));
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
%                     [road, systemState, growth_energy, vdir,uselessRuns,moved] = Walking(systemState, jj, prob_lateral, prob_forward, prob_lateralp, prob_forwardp, walk_ECost,growth_energy,m, Transp_ECost, road, x, y, roadlength, tunneltip, elow_limit, prob_turn,tunnelsize,uselessRuns,countSpot,moved);
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
                % image displays only positive values: +2 rescales
                %                 figure(1);
                %                 image(road+2);
                %                 colormap([0.25 0.25 0.25; 0.8 0.8 0; 1 1 1; 1 0.5 0]);
                %                 drawnow;
                
                figure(2);
                clf;
                set(gca,'box','on');
                
                newR = road+2;
                for i = 1:numantsALL
                    id = find(systemState.globalindx==i);
                    %                     if(systemState.digpause(id)>1)
                    %                         newR(systemState.y(id),...
                    %                             systemState.x(id))=1;
                    %                     else
                    newR(systemState.y(id),...
                        systemState.x(id))=i+4;
                    %                     end
                    
                end
                %                 cmap = lines(numantsALL+4);
                image(newR);
                
                text(10,1,['t = ',num2str(kk/2/60,'%10.2g'),' min'],'fontsize',20);
                colormap(cmap);
                colorbar;
                drawnow;
                %                 pause;
                
                %                 pause;
                %pause;
                %
                
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
            
            k1=find(abs(road(1,(roadlength-tunneltip+countSize):(roadlength-1)))==2);
            k2=find(abs(road(2,(roadlength-tunneltip+countSize):(roadlength-1)))==2);
            
            density(kk,:)= [length(k1),length(k2)];
            if tunneltip==1
                break;
            end
            
            %             probabil(systemState.globalindx(:))= systemState.prob(:);
            %             atFace(systemState.globalindx(:))=atFace(:);
            
            
            restInds= find(systemState.restflag-prevResting==-1); %resting=1 to walking=0 0-1=-1
            resting(restInds)=resting(restInds)+1;
            prevResting=systemState.restflag(:);
            
        end
        
        %after diff==-1 means pellet deposited
        
        pellet =[zeros(1,numants); cumsum(diff(pellet)==-1)];
        
        % set energy
        growth_energy=growth_energy(m,:);
        
%         groupEnergy(:,:,run)=growth_energy(1:size(growth_energy),:); % total energy expenditure per unit of tunnel growth: 1 -tunneltip; 2 -energy expenditure; 3 - number of steps kk
        indEnergy(:,:,run)= energy; % energy of individual ants per 1 step of cycle
        tunLength(:,:,run)=tunnel_length; % tunnel growth
        %         markMatr(:,:)=markmatr;
        %         markMatrN0(:,:,run)=markmatrN0;
        tunTime(systemState.globalindx(:))=tunTime(:);
        
    end
    
    ResEnergy(z).numants=numants;
    ResEnergy(z).groupEnergy=groupEnergy;
        ResEnergy(z).groupEnergy=growth_energy;
    ResEnergy(z).indEnergy=indEnergy;
    ResEnergy(z).tunLength=tunLength;
    %         ResEnergy(schRes).markMatr=markMatr;   %%? next line it changes?
    ResEnergy(z).density=density;
    ResEnergy(z).markMatr=markmatrN0;
    ResEnergy(z).pall1=all1;
    ResEnergy(z).coeffP=coeffP;
    ResEnergy(z).penalty=penalty;
    ResEnergy(z).uselessRuns= uselessRuns;
    ResEnergy(z).numants=numants;
    ResEnergy(z).pellets=pellet;
    ResEnergy(z).prob=probs;
    ResEnergy(z).atFace= atFace;
    ResEnergy(z).resting=resting;
    ResEnergy(z).tunTime = tunTime;
    ResEnergy(z).flow = flow;
    ResEnergy(z).prob_turn = prob_turn;
    ResEnergy(z).pause2dig=pause2dig;
    ResEnergy(z).recharge_steps = recharge_steps;
    ResEnergy(z).occupied = occupied;
    ResEnergy(z).equalDis     = all(probs(1,1)==probs(1,:));
    ResEnergy(z).pell2grow = pellet2grow;
    ResEnergy(z).infEnergy = energyMult>1;
    ResEnergy(z).equalDis     = all(probs(1,1)==probs(1,:));
    %energyMult = 1 if non-inf energy!
%     schRes=z+1;
    clearvars groupEnergy indEnergy tunLength markMatr markMatrN0
end
%trim last line in group Energy
% save([pwd,'\ResEnergyP0ESame.mat'],'ResEnergy', '-v7.3')
% if(nargin==0)
%     lorenzcurve(ResEnergy.prob,ResEnergy.pellets);
% end
% if nargin<6
toc
% end