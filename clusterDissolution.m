%This program finds out cluster dissolution times for clusters of various
%dimensions %originally written by Bahnisikha Dutta
clc;
clearvars -except road res r

% TestImages=load('TestImages.mat');%Load CA Frames
% Cluster_Info=clusterInfo(TestImages.TestImages); %Refer to function ClusterInformation


Cluster_Info=clusterInfo(road); %Refer to function ClusterInformation

%% Initialization
frame_no=length(Cluster_Info);
max_clusters=numel(Cluster_Info)/length(Cluster_Info);
com_x=zeros(max_clusters,frame_no);
com_y=zeros(max_clusters,frame_no);
size=zeros(max_clusters,frame_no);
temp_y=zeros(max_clusters,frame_no-1);
temp_x=zeros(max_clusters,frame_no-1);
ones_y=zeros(max_clusters,frame_no-1);
ones_x=zeros(max_clusters,frame_no-1);
dissolution_time=zeros(1,max_clusters);

%% Retrieve size and centre of mass of each cluster
for i=1:1:frame_no
    for j=1:1:max_clusters
        if (isempty(Cluster_Info(j,i).size)==1)
             com_x(j,i)=0;
             com_y(j,i)=0;
             size(j,i)=0;
        else
             com_x(j,i)=Cluster_Info(j,i).medianx;
             com_y(j,i)=Cluster_Info(j,i).mediany;
             size(j,i)=Cluster_Info(j,i).size;
        end 
    end
end

%% Check for change in centre of mass for each cluster
for i=1:1:frame_no-1

        temp_x1=com_x(:,i);
        temp_x2=com_x(:,i+1);
        temp_y1=com_y(:,i);
        temp_y2=com_y(:,i+1);
        temp_x(:,i)=temp_x2-temp_x1;
        temp_y(:,i)=temp_y2-temp_y1;

end

        [mask_xi,mask_xj]=find(com_x==0);
        [mask_yi,mask_yj]=find(com_y==0);
        for k=1:1:numel(mask_yi)
        ones_masky(mask_yi(k),mask_yj(k))=100;
        end
        for l=1:1:numel(mask_xi)
        ones_maskx(mask_xi(l),mask_xj(l))=100;
        end
        
         temp_y=temp_y+ones_masky(:,2:end);
         temp_x=temp_x+ones_maskx(:,2:end);
        [temp_yi,temp_yj]=find(temp_y==0);
        [temp_xi,temp_xj]=find(temp_x==0);
        
%% Create a mask of 1's where centre of mass did not change        
        for k=1:1:numel(temp_yi)
        ones_y(temp_yi(k),temp_yj(k))=1;
        end
        for l=1:1:numel(temp_xi)
        ones_x(temp_xi(l),temp_xj(l))=1;
        end
        ones=ones_x.*ones_y; 
%% Calculate dissolution time for each cluster and link it to cluster size         
  for m=1:1:max_clusters
      ones_ind=find(ones(m,:)==1);
      diff_ones=diff(ones_ind);
      if(find(diff_ones>1)>0)
      cluster_no=find(diff_ones>1);
      cluster_size=diff(cluster_no);
      dissolution_time(1)=cluster_no(1);
      dissolution_time(2:(numel(cluster_size)+1))=cluster_size;
      dissolution_time(end+1)=numel(ones_ind)-cluster_no(end);
      cluster_ind(1)=ones_ind(1);
      for j=2:1:(numel(cluster_no)+1)
      cluster_ind(j)=cluster_ind(j-1)+diff_ones(cluster_no);
      end
      
      for i=1:1:(numel(cluster_ind))
      rel_cluster_size(i)=size(m,cluster_ind(i));
      end
      
      else
          dissolution_time=numel(diff_ones)+1;
          rel_cluster_size=size(m,ones_ind(1));
      end
      dissolution_time_cell{m}=dissolution_time;
      cluster_size_cell{m}=rel_cluster_size;
      clear dissolution_time rel_cluster_size
  end
  
  %%Plot the cluster sizes against resolution time
  time=cell2mat(dissolution_time_cell);
  c_size=cell2mat(cluster_size_cell);
  figure()
  plot(c_size,time,'k*');
  axis equal
  xlabel('Cluster Sizes');
  ylabel('Dissolution Time');
        