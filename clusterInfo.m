function [ cluster_data ] = Cluster_Information( TestImages )
%Author:Bahnisikha Dutta
%Date:08/14/2017
%Cluster_Info function returns the information about clusters present in every
%CA frame 
%Input arguement is a 3d matrix of the CA frames
% Returned Structure contains cluster informations including size of
% clusters, coordinates of clusters and median of clusters
for j=1:1:size(TestImages,3)
test=TestImages(:,:,j);%frame to be analysed
[onesx,onesy]=find(test==1);%this finds out indices for all the ants
%possible exclusions: a)just one ant present b)diagonal ants c)just one
%cluster of multiple ants d)multiple one ant clusters e)2 or more ants
%sharing the same y coordinate(or following a line)-considering this still
%a cluster; will determine if its a clog or not while postprocessing center
%of mass
if ((numel(onesx)>1)||(numel(onesy)>1)) %this checks if there are more than one ant-case (a)
    
diff_y=diff(onesy);% first differential to find out the breaks in coordinates to demarcate clusters

if(find(diff_y>1)>0) %checks if there are more than one clusters- case (c)
cluster_no=find(diff_y>1);
cluster_size=diff(cluster_no);
clusters(1)=cluster_no(1);
clusters(2:(size(cluster_size)+1))=cluster_size;
clusters(end+1)=size(onesy,1)-cluster_no(end);
% cluster_no(end+1)=size(onesy,1);
k=1;
for i=1:size(clusters,2)
    if (clusters(i)>1)
        %this takes care of diagonal ant cases-case (b)-this needs to be
        %modified if tunnel width is greater than 2
        if (clusters(i)==2)
            check_diag=onesx(k:k-1+clusters(i));
            diff_diag=diff(check_diag);
            if(diff_diag==0)
            clusters_structure(i,j).coordinatex=onesx(k:k-1+clusters(i));
            clusters_structure(i,j).coordinatey=onesy(k:k-1+clusters(i));
            clusters_structure(i,j).size=clusters(i);
            clusters_structure(i,j).medianx=median(onesx(k:k-1+clusters(i)));
            clusters_structure(i,j).mediany=median(onesy(k:k-1+clusters(i)));
            end
        else 
            clusters_structure(i,j).coordinatex=onesx(k:k-1+clusters(i));
            clusters_structure(i,j).coordinatey=onesy(k:k-1+clusters(i));
            clusters_structure(i,j).size=clusters(i);
            clusters_structure(i,j).medianx=median(onesx(k:k-1+clusters(i)));
            clusters_structure(i,j).mediany=median(onesy(k:k-1+clusters(i)));
        end
    end
    k=k+clusters(i);
end


else
    clusters_structure(1,j).coordinatex=onesx;
    clusters_structure(1,j).coordinatey=onesy;
    clusters_structure(1,j).size=size(onesx,1);
    clusters_structure(1,j).medianx=median(onesx);
    clusters_structure(1,j).mediany=median(onesy);
end
end
clear cluster_size cluster_no clusters onesx onesy %test check_diag diff_diag diff_y
end
cluster_data=clusters_structure;

end

