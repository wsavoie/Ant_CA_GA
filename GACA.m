for iii=([0 2])
    clearvars -except iii
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
    if(isunix)
        outputFolder='/home/ws/Documents/Ant_CA_GA';
    end
    ObjectiveFunction = @GA_CA_code;
    
    
    %%%%%rech and ptt%%%%%
    % nvars = 32;
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%regular%%%%%%%%
    nvars = 30;
    %25 30 35 40 45 50 55 60 65 70
    %%%%%%%%%%%%%%%%%%%%%%
    totgens = 50;
    popsize = 5*nvars;
    tuntip=10;
    prob2turn = .3;
    rechargeSteps = 600;
    numIts=432*4;
    energyMult=0;
    TW=2;
    mutType=1; %0=adapt 1=uniform
    mutRate=.4;
    LB = zeros(1,nvars);
    UB = ones(1,nvars);
    bestofgen= cell(1,totgens);
    initSeed=iii; %0=equal 1=rand 2=unequal
    
    
    crossOverFrac = 0; %no crossover
    % options = optimoptions(@ga,'MutationFcn',{@mutationgaussian,.2,.5});
    if mutType==0
        options = optimoptions(@ga,'MutationFcn', {@mutationadaptfeasible});
    elseif mutType==1
        options = optimoptions(@ga,'MutationFcn', {@mutationuniform, mutRate});
    end
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
    
    if initSeed==0
        %%%%%%%%EQUAL%%%%%%%%%%%%%%%
        mypop= ones(popsize,nvars);
        options = optimoptions(options,'InitialPopulation', mypop);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif initSeed==1
        %%%%%%%%RANDOM%%%%%%%%%%%%%%%
        mypop= rand(popsize,nvars);
        options = optimoptions(options,'InitialPopulation', mypop);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif initSeed==2
        %%%%%UNEQUAL%%%%%%%
        mypop= zeros(popsize,nvars);
        mypop(:,nvars)=1;
        options = optimoptions(options,'InitialPopulation', mypop);
        %%%%%%%%%%%%%%%%%%%
    end
    
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
    fname=['N=',num2str(nvars),'_','tw=',num2str(TW),'_',char(datetime('now','Format','yyyy-MM-dd-HH-mm'))];
    clear myCluster
    save(fullfile(outputFolder,'results',[fname,'.mat']));
    paramsFilePath=fullfile(outputFolder,'results',[fname,'.txt']);
    fileID = fopen(paramsFilePath,'w');
    if(isunix)
        outputType='%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n';
    else
        outputType='%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n';
    end
    fprintf(fileID,outputType,...
        ['nvars = ',num2str(nvars),';'],...
        ['totgens = ',num2str(totgens),';'],...
        ['popsize = ',num2str(popsize),';'],...
        ['tuntip = ',num2str(tuntip),';'],...
        ['prob2turn = ',num2str(prob2turn),';'],...
        ['rechargeSteps = ',num2str(rechargeSteps),';'],...
        ['numIts = ',num2str(numIts),';'],...
        ['energyMult = ',num2str(energyMult),';'],...
        ['TW = ',num2str(TW),';'],...
        ['MutType = ',num2str(mutType),';'],...
        ['MutRate = ',num2str(mutRate),';'],...
        ['initSeed= ',num2str(initSeed),';']);
    
    fclose(fileID);
    close all;
    
end