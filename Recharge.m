function [road, systemState] = Recharge(road, systemState, recharge, etop_limit, x, y, jj,coeffP)
%RECHARGE Summary of this function goes here
%-2 down
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

