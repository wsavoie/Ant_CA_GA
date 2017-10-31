%This program finds out cluster dissolution times for clusters of various;
%case: identity of ants not preserved
%dimensions
clc;
clear all;


Ca_equal_mat=[0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]; %uncomment while bulk processing only

% figure()
hold on
% for z=1:8:numel(Ca_equal_mat)%uncomment while bulk processing only
%  T='type_1_R_0.01.mat';
T='type_0_R_0.8.mat';
 fold='D:\Projects\Ant_CA_GA\results\longRuns 50 gens recharge .4 mut\finEng_24h\reversal data 432 its'
%  fold='D:\Projects\Ant_CA_GA\results\older 10-17\bahni workload old';
% T=strcat('type_1_R_',num2str(Ca_equal_mat(z)), '.mat');%uncomment while bulk processing only
TestImages=load(fullfile(fold,T));%uncomment while bulk processing only
% TestImages=load('roadFull.mat');%Load CA Frames,comment while bulk processing
Cluster_Info=Cluster_Information(TestImages.roadFull);%uncomment while bulk processing only
% Cluster_Info=Cluster_Information(TestImages.road); %Refer to function,comment while bulk processing
% ClusterInformation
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
      if (numel(ones_ind)==0)
      dissolution_time_cell{m}=0;
      cluster_size_cell{m}=0;
      else
      diff_ones=diff(ones_ind);
     
      if(find(diff_ones>1)>0)
            cluster_no=find(diff_ones>1);
            cluster_size=diff(cluster_no);
            dissolution_time(1)=cluster_no(1);
            dissolution_time(2:(numel(cluster_size)+1))=cluster_size;
            dissolution_time(end+1)=numel(ones_ind)-cluster_no(end);
            cluster_ind(1)=ones_ind(1);
            for j=2:1:(numel(cluster_no)+1)
                cluster_ind(j)=cluster_ind(j-1)+diff_ones(cluster_no(j-1));
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
      end
      clear dissolution_time rel_cluster_size cluster_ind
  end
  
  %%Plot the cluster sizes against resolution time
  time_temp=cell2mat(dissolution_time_cell);
  c_size_temp=cell2mat(cluster_size_cell);
  c_size=c_size_temp(find(c_size_temp>0));
  time=time_temp(find(c_size_temp>0));
%   %% comment while bulk processing
%   figure(1)
%   plot(c_size,time,'k*');
%   axis equal
%   xlabel('Cluster Sizes');
%   ylabel('Dissolution Time');
  
  %%
  c_size_un=unique(c_size); 
  for i=1:1:numel(c_size_un)
      ind=find(c_size==c_size_un(i));
      freq(i)=numel(ind);
      time_ind=time(ind);
      mean_ind(i)=mean(time_ind);
      std_ind(i)=std(time_ind);
  end
  
%     freq=freq/max(freq);
%     mean_ind=mean_ind/max(mean_ind);
%     std_ind=std_ind/max(std_ind);
    
% %%   comment while bulk processing 
%     figure()
%     
%    [hAx,hLine1,hLine2]=plotyy(c_size_un,freq,c_size_un,mean_ind);
%     xlim(hAx(1),[min(c_size_un)-0.5 max(c_size_un)+0.5]);
%    hLine1.LineStyle = '--';
%    hLine2.LineStyle = ':';
%     hold on
%     bar(c_size_un,freq);
% %     ylabel('Normalized Frequency');
%     hold(hAx(2), 'on');
%     ylim(hAx(2),[-1 3]);
%     errorbar(hAx(2), c_size_un,mean_ind,std_ind);
% %     ylabel('Normalized Dissolution Time');
%     xlabel('Cluster Sizes');
%     ylabel(hAx(1),'Normalized Frequency','FontSize',12, 'FontWeight','Bold'); % left y-axis 
%     ylabel(hAx(2),'Normalized Dissolution Time','FontSize',12, 'FontWeight','Bold')
%%   %uncomment while bulk processing only
%     figure()
    plot(c_size_un,freq/sum(freq),'LineWidth',2);
    xlabel('Cluster Sizes');
    ylabel('Frequency'); 
%     set(gca,'FontSize',18, 'FontWeight','Bold', 'LineWidth',2);
figText(gcf,16);
    hold on;
    clearvars -except Ca_equal_mat ;%uncomment while bulk processing only
% end%uncomment while bulk processing only
%     legend('unequal R=0.3', 'equal R=0.3');

% legend('equal R=0.01', 'equal R=0.1', 'equal R=0.2', 'equal R=0.3', 'equal R=0.4', 'equal R=0.5', 'equal R=0.6', 'equal R=0.7', 'equal R=0.8'); %uncomment while bulk processing only
% legend('equal R=0.01', 'equal R=0.1', 'equal R=0.2', 'equal R=0.3', 'equal R=0.4', 'equal R=0.5', 'equal R=0.6', 'equal R=0.7', 'equal R=0.8','unequal R=0.01', 'unequal R=0.1', 'unequal R=0.2', 'unequal R=0.3', 'unequal R=0.4', 'unequal R=0.5', 'unequal R=0.6', 'unequal R=0.7', 'unequal R=0.8');
% legend('unequal R=0.01', 'unequal R=0.1', 'unequal R=0.2', 'unequal R=0.3', 'unequal R=0.4', 'unequal R=0.5', 'unequal R=0.6', 'unequal R=0.7', 'unequal R=0.8'); %uncomment while bulk processing only


%%Plots for excavation efficiency v/s reversal probabilities
load('reversalVals.mat');
equal_ind=find(datAll(:,1)==0);
equal_revprob=datAll(equal_ind,4);
unequal_ind=find(datAll(:,1)==1);
unequal_revprob=datAll(unequal_ind,4);
rev_prob=[0.0100000000000000;0.0200000000000000;0.0300000000000000;0.0400000000000000;0.0500000000000000;0.0600000000000000;0.0700000000000000;0.0800000000000000;0.0900000000000000;0.100000000000000;0.125000000000000;0.150000000000000;0.175000000000000;0.200000000000000;0.225000000000000;0.250000000000000;0.275000000000000;0.300000000000000;0.325000000000000;0.350000000000000;0.375000000000000;0.400000000000000;0.425000000000000;0.450000000000000;0.475000000000000;0.500000000000000;0.525000000000000;0.550000000000000;0.575000000000000;0.600000000000000;0.625000000000000;0.650000000000000;0.675000000000000;0.700000000000000;0.725000000000000;0.750000000000000;0.775000000000000;0.800000000000000;0.825000000000000;0.850000000000000;0.875000000000000;0.900000000000000;0.925000000000000;0.950000000000000;0.975000000000000;1];

figure()

for i=1:1:numel(rev_prob)
plot(rev_prob(i),equal_revprob(i),'o','LineWidth',2,'MarkerSize',10);
hold on;
end

figure()

for i=1:1:numel(rev_prob)
plot(rev_prob(i),unequal_revprob(i),'o','LineWidth',2,'MarkerSize',10);
hold on;
end

% legend('R=0.01', 'R=0.1', 'R=0.2', 'R=0.3', 'R=0.4', 'R=0.5', 'R=0.6', 'R=0.7', 'R=0.8'); %uncomment while bulk processing only
xlabel('Reversal probabilities','FontSize',36, 'FontWeight','Bold');
ylabel('Excavation Efficiency(cm)','FontSize',36, 'FontWeight','Bold');
set(gca,'FontSize',36, 'FontWeight','Bold','LineWidth',2);
title('Equal');
title('Lorenz/Unequal');