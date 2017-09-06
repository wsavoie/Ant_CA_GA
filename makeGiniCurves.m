%ugly way of finding random gini coefficients

gMin=.1;
gRes=.1;
thresh = .001;
totcc=1; %total iterations through loop
cc=1;    %iterations through ginis
ants=30;

giniX=gMin:gRes:.9;
gX=giniX; %stores all remaining gini options

giniN=length(giniX);
g=zeros(giniN,1);
giniY=zeros(giniN,ants);

h=waitbar(0,[num2str(0),'/',num2str(giniN),' total iterations=',num2str(totcc)]);

while ~isempty(gX)
    totcc=totcc+1;
    %generate random probability
    p=rand(ants,1);
    p=p/sum(p);
    [G,~]=Gini(p);
    
    v=gX(G<(gX(G>(gX-thresh))+thresh));
    if(v)
        gX(gX==v)=[];
        ind=find(giniX==v);
        waitbar(cc/giniN,h,[num2str(cc),'/',num2str(giniN),' total iterations=',num2str(totcc)]);
        giniY(ind,:)=p;
        cc=cc+1;
        for i=1:giniN 
            g(i)=Gini(giniY(i,:)); 
        end
    end
    if(mod(totcc,100000)==0)
        pts('total iterations:', totcc);
    end
end
closeWaitbar
%it likely will get stuck running for a long time, so make last few using
%below code:
x=17;p(x:30)=[rand(30-x+1,1)]; p=p/sum(p); Gini(p);

% for(i=1:giniN) GG(i)=Gini(gGiniY(i,:)); end
% for(i=1:giniN) g(i)=Gini(giniY(i,:)); end