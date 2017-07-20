% function [systemState,road,tunneltip,growth_energy,m,kk,pelletcount,pellet2grow] = Digging(systemState,jj,pause2dig,exc_ECost,walk_ECost,growth_energy,m, road, x, y, roadlength, tunneltip, kk, pelletcount, pellet2grow)
function [systemState,road,tunneltip,growth_energy,m,kk,pelletcount,pellet2grow] = Digging(systemState,jj,pause2dig,exc_ECost,~,growth_energy,m, road, x, y, roadlength, tunneltip, kk, pelletcount, pellet2grow)

%DIGGING Summary of this function goes here
if systemState.digpause(jj,1)==pause2dig
    systemState.digpause(jj,1)=0;
    systemState.pellet(jj,1)=1;
    systemState.energy(jj,1)=systemState.energy(jj,1)-exc_ECost;
    growth_energy(m,2)=growth_energy(m,2)+exc_ECost;
    pelletcount=pelletcount+1;
    %remove grain at the pausetodig step
    if pelletcount==pellet2grow
        pelletcount=0;
        xx= find(road(:,roadlength-tunneltip-1)==0,1);
        if ~isempty(xx) %there exists a space along tunnel tip with sand
            if road(y(jj),roadlength-tunneltip-1)==0
                road(y(jj),roadlength-tunneltip-1)=1;
            else
                road(xx,roadlength-tunneltip-1)=1;
            end
        end
        xx= find(road(:,roadlength-tunneltip-1)==0,1);
        if isempty(xx)
            %    if isempty(find(road(:,roadlength-tunneltip-1)==0))
            tunneltip=tunneltip+1;
            %HERE ROAD INCREASES IN LENGTH -> growth_energy(m+1,1)=0;
            growth_energy(m+1,2)=growth_energy(m,2);
            %disp(tunneltip);
            growth_energy(m,1)=tunneltip;
            growth_energy(m,3)=kk; % number of iteration
            m=m+1;
        end
        
    end
    
    %change direction@ tip
    road(y(jj), x(jj))=-road(y(jj),x(jj));
else
    % pause for digging
    systemState.digpause(jj,1)=systemState.digpause(jj,1)+ 1;% DOUBLE CHECK THIS!! +1
    systemState.energy(jj,1)= systemState.energy(jj,1)-exc_ECost;
    growth_energy(m,2)=growth_energy(m,2)+exc_ECost;
end


end

