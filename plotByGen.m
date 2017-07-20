maxpells=zeros(length(gens),1);
ind=zeros(length(gens(1).out),1);
g
for i=1:length(gens)
    for j = 1:length(gens(i).out)
         if(maxpells(i) < sum(gens(i).out(j).pellets))
            maxpells(i) = sum(gens(i).out(j).pellets);
            ind(i) = j;
         end
    end
    len = length(gens(i).out(ind(i)).pellets/sum(gens(i).out(ind(i)).pellets));
    dat =gens(i).out(ind(i)).pellets/sum(gens(i).out(ind(i)).pellets);
   [g(i)]=lorenzcurve(ones(len,1)/len,dat);
   close 1
end

figure(20)
hold on;
title('gini coefficient')
xlabel('generations')
plot(1:length(gens),g);
figure(21)
hold on;
title('pellets excavated')
xlabel('generations')
plot(1:length(gens),maxpells(1:length(gens)));
