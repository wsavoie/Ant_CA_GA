function [G,xy]=Gini(y,varargin) 


if(~iscolumn(y))
    y=y';
end
x = ones(length(y),1)/length(y);

y = y/ sum(y);
xy = [x, y, y ./ x];
% Sort with respect to Wealth Per Capita
xy = sortrows(xy, 3);
xy(:, 3) = [];
% Cumulative p & w
xy = cumsum(xy);
xy = [zeros(1, 2); xy];
A=polyarea([xy(:,1);xy(1,1)],[xy(:,2);xy(1,2)]); 
B=polyarea([xy(:,1);xy(end,1);xy(1,1)],[xy(:,2);xy(1,2);xy(1,2)]);

G=A/(A+B);
if(nargin>1)
    hold on;
    fill([xy(:,1);xy(1,1)],[xy(:,2);xy(1,2)],'y','FaceColor','y','EdgeColor','y'); %area a;
    fill([xy(:,1);xy(end,1);xy(1,1)],[xy(:,2);xy(1,2);xy(1,2)],'b','FaceColor','b');
    plot(xy(:,1),xy(:,2),'k','LineWidth',4);
    axis([0 1 0 1])
    figText(gcf,15);
    legend({'Area A','Area B',horzcat('Lorenz Curve G=',num2str(G))})

end

end