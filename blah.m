gy(:,1),gy(:,2)
x1=gy(:,1);
y1=gy(:,2);
x2=[0,1];
y2=[1,0];
px2 = interp1(x2,y2,x1,'linear');
py2 = interp1(y2,y2,x1,'linear');
% xx = xx(max(x1(1),x2(1)) <= xx & min(x1(end),x2(end)) >= xx);
func = @(x)ppval(px2,x)-ppval(py2,x);
xb = xx([true; diff(func(xx) > 0) ~= 0]);
i1 = hankel(1:2,2:numel(xb));
xout = arrayfun(@(z)fzero(func, xb(i1(:,z))), (1:size(i1,2))' )