function [road, systemState, growth_energy, rev,uselessRuns,moved] = Walking(systemState, jj, prob_lateral, prob_forward, prob_lateralp, prob_forwardp, walk_ECost,growth_energy,m, Transp_ECost, road, x, y, roadlength, tunneltip, elow_limit, prob_turn,tunnelsize,uselessRuns,countSpot,moved)
%WALKING Summary of this function goes here
xorig= x(jj);
yorig= y(jj);
rev=0;
%decide if ant carries a pellet or not

if  systemState.pellet(jj,1)==0
    lateralmove = rand(1) < prob_lateral;
    forwardmove = rand(1) < prob_forward;
    lateralmove = lateralmove*(2*(rand(1) < 0.5)-1);% lateral move can be -1,1
    systemState.energy(jj,1)=systemState.energy(jj,1)-walk_ECost;
    growth_energy(m,2)=growth_energy(m,2)+walk_ECost;
else
    lateralmove = rand < prob_lateralp;
    forwardmove = rand < prob_forwardp;
    lateralmove = lateralmove*(2*(rand(1) < 0.5)-1);% lateral move can be -1,1
    systemState.energy(jj,1)=systemState.energy(jj,1)-Transp_ECost;
    growth_energy(m,2)=growth_energy(m,2)+Transp_ECost;
end

%ant randomly changes the lane
%boundary conditions
if(lateralmove && y(jj) == 1)
    lateralmove = double(1);
end

if(lateralmove && y(jj) == tunnelsize)
    lateralmove = double(-1);
end

% figure out direction of motion +1 up /-1 down
vdir = road(y(jj), x(jj))/abs(road(y(jj), x(jj)));
%forvardmove=+1, 0 - next step - stop/move forward

%predict next ant position
xahead = x(jj) + forwardmove*vdir;
yahead = y(jj) + lateralmove;

% HERE WE CHECK IF ANTS HAVE ENOUGH ENERGY FOR NEXT TRIP & DECIDE WHAT ANT
% DOES AT THE EXIT
if(x(jj) == roadlength && vdir > 0) % +1 - up
    if systemState.pellet(jj,1)==0
        uselessRuns=uselessRuns+1;
    end
    
    %do nothing for 1 step and change direction @ exit
    if systemState.energy(jj,1)<elow_limit
        %Flag to stop digging
        systemState.restflag(jj,1)=1;
        xahead = roadlength;
        systemState.pellet(jj,1)=0;
        % NO ENERGY-ANT DROPS PELLET AND GOES TO RECHARGE
        % ant vanishes until recharged
        road(y(jj), x(jj))=1;
    else
        % ant left the tunnel but did not reach low
        % energy state
        systemState.restflag(jj,1)=1;%rest flag tells that ant is out of the tunnel
        xahead = roadlength;
        systemState.pellet(jj,1)=1;% ant still has a pellet and will drop it after X-deposition steps
        %ant exits to drop a pellet EVEN IF ANTS HAS NO PELLET WE GAVE IT A
        %PELLET HERE SO IT WOULD SPEND 20 STEPS AT THE TOP!!! AS A RESULT
        %EVERY ANT PAUSES BEFORE SHE GOES BACK!
        road(y(jj), x(jj))=1;
    end
    
    
end

%random turn back when descending ant meets ascending one % ANTS WHO GOES
%UP CAN ALSO CHANGE DIRECTION-?

if abs(road(yahead, xahead)) == 2 && xahead>=roadlength-tunneltip && vdir<0 && xahead~=x(jj)&& systemState.restflag(jj,1)~=1
    %disp('here')
    if rand(1)<prob_turn
        road(y(jj), x(jj)) = -road(y(jj), x(jj)); % random turn when jammed
        rev=1;
    else
        rev=0;
    end
end

% now check if occupied IF EMPTY - MOVE AHEAD


if (road(yahead, xahead) == 1 && xahead>=roadlength-tunneltip)
    road(yahead, xahead) = road(y(jj), x(jj));
    road(y(jj), x(jj)) = 1;
    systemState.x(jj,1)=xahead;
    systemState.y(jj,1)=yahead;
    
elseif abs(road(yahead, xahead)) == 2
    %do not move if there is someone in front of you
    systemState.x(jj,1)=x(jj);
    systemState.y(jj,1)=y(jj);
    
end

if(xorig==countSpot)
%     if(xorig ~= systemState.x(jj)||yorig~=systemState.y(jj)) %orig
    if((xorig ~= systemState.x(jj)||yorig~=systemState.y(jj))&& systemState.pellet(jj)==1)
        moved(yorig)=1;
    else
        moved(yorig)=0;
    end
end