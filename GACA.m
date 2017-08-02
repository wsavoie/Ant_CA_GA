
% ObjectiveFunction = @simple_fitness;
% nvars = 4;    % Number of variables
% LB = [0 0 0 0];   % Lower bound
% UB = [30 30 30 30];  % Upper bound
% ConstraintFunction = @simple_constraint;
%
% options = optimoptions('PlotFcn');
% options = optimoptions(options,'PlotFcn',{@gaplotbestf,@gaplotmaxconstr}, ...
%     'Display','iter');
% % Next we run the GA solver.
% [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ...
%     ConstraintFunction,options)\
tic
outputFolder='D:\Projects\Ant_CA_GA';

ObjectiveFunction = @GA_CA_code;


%%%%%rech and ptt%%%%%
% nvars = 32;
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%regular%%%%%%%%
nvars = 100;
%%%%%%%%%%%%%%%%%%%%%%
totgens = 20;
popsize = 3;
tuntip=10;
prob2turn = .3;
rechargeSteps = 10;
numIts=432;
energyMult=1;
TW=5;

LB = zeros(1,nvars);
UB = ones(1,nvars);
bestofgen= cell(1,totgens);



crossOverFrac = 0; %no crossover
% options = optimoptions(@ga,'MutationFcn',{@mutationuniform,mutRate});
options = optimoptions(@ga,'MutationFcn', {@mutationadaptfeasible});

%%%%%%%%%%%%%%%%%%%%normal%%%%%%%%%%%%%%%%%%%%
options = optimoptions(options,'PlotFcn',{{@outputFunc,bestofgen},...
    {@gaplotlorenzcurve,bestofgen,numIts,rechargeSteps,prob2turn,tuntip,TW,energyMult}}, ...
    'Display','iter','MaxGenerations',totgens,'CrossoverFraction',crossOverFrac);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%rech and ptt%%%%%%%%%%%%%%%%%%%%
% LB(31:32) = [100 .3];
% UB(31:32)=[800 .8];
% options = optimoptions(options,'PlotFcn',{{@outputFuncTest,bestofgen},...
%     {@gaplotlorenzcurve,bestofgen,numIts}}, ...
%     'Display','iter','MaxGenerations',totgens);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%EQUAL%%%%%%%%%%%%%%%
% mypop= ones(popsize,nvars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%UNEQUAL%%%%%%%
% mypop= zeros(popsize,nvars);
% mypop(:,nvars)=1;
%%%%%%%%%%%%%%%%%%%
% options = optimoptions(options,'InitialPopulation', mypop);

options = optimoptions(options,'PopulationSize', popsize);

myCluster=parcluster('local');
myCluster.NumWorkers=7;
saveProfile(myCluster);
options = optimoptions(options,'UseParallel', true);
% proba,tt,rech,ptt,tuntip,tw,energyMult
%%%NORMAL%%%
[x,fval,exitflag,output,population,scores]=ga({ObjectiveFunction,numIts,rechargeSteps,prob2turn,tuntip,TW,energyMult},nvars,[],[],[],[],LB,UB,[],options);
%%%%RECH AND PTT%%%
% [x,fval,exitflag,output,population,scores]=ga({ObjectiveFunction,numIts},nvars,[],[],[],[],LB,UB,[],options);
toc

%%
bestofgen = bestofgen(~cellfun(@isempty,bestofgen));
bestofgenOUT=cell(1,length(bestofgen));
% y=bestofgen{end};
clear res;
parfor i=1:length(bestofgen)
    y=bestofgen{i};
    
    % res(i) =CA_FunctionsWill(y,length(y),numIts,2,0);  %probs,numants,numits*10000,width,infEnergy
    % f=sum(res(i).markMatr(:,2:end));
    % res =CA_FunctionsWill(y,length(y),numIts,tw,emult,1,rechargeSteps,prob2turn);  %probs,numants,numits*10000,width,infEnergy
    res=CA_FunctionsWill(y,length(y),numIts,TW,energyMult,1,rechargeSteps,prob2turn,tuntip);
    f=sum(res.markMatr(:,2:end));
    bestofgenOUT{i}=f;
end

% save(['A:\','PAPERFIG_RAND_R=',num2str(rechargeSteps),'P=',num2str(prob2turn),'.mat']);
% date
save(fullfile(outputFolder,'results',['50ants_',char(datetime('now','Format','yyyy-MM-dd-HH-mm')),'.mat']));