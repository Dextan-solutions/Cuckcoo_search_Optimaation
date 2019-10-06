clc;clear;close all;

%% User Input
DistLoadFlowSolution=powerflow;
% Function Input
User.Function='ObjfuncPoweRnLSF';
User.NumbVar=4;
% User.Lb=[-20,-10,-5,-1];
% User.Ub=[20,10,5,1];

User.MaxIter=100;
User.Lb=[1000 1];
User.Ub=[5000 33];

% User.Lb=[0 0 0];
% User.Ub=[2.0  2  2];
% Cuckoo Input
User.NumNest=25;
%Levy Flight Input
User.beta=3/2;
pa=0.25;

    %% Initializing the Cuckoo Algorithm
SampleNest.DGnLoc=[];
SampleNest.PLosVolt=[];
Nest=repmat(SampleNest,User.NumNest,1);
for i=1:User.NumNest
%     Nest(i).Position=round(User.Lb+(User.Ub-User.Lb).*rand(size(User.Lb)));
    Nest(i).DGnLoc=round(User.Lb+(User.Ub-User.Lb).*rand(size(User.Lb)));
    
       DistLoadFlowDGSolution=powerflowDG(Nest(i).DGnLoc(1,1),Nest(i).DGnLoc(1,2));
    
       Nest(i).CostPLos=[DistLoadFlowDGSolution.PtLosskW];
       Nest(i).CostPbrLos=[DistLoadFlowDGSolution.Pbrloss];
       Nest(i).CostVact=[DistLoadFlowDGSolution.Vactual];
       Nest(i).CostVolt=[DistLoadFlowDGSolution.VmagPU];
       Nest(i).CostVSI=[DistLoadFlowDGSolution.VSI];
       Nest(i).CostMinVolt=[DistLoadFlowDGSolution.minVSI];
       Nest(i).CostLSF=[DistLoadFlowDGSolution.LSF];
       Nest(i).CostVangle=[DistLoadFlowDGSolution.Vangle];
       Nest(i).CostQtLos=[DistLoadFlowDGSolution.QtLosskVAr];
       Nest(i).CostQbrLos=[DistLoadFlowDGSolution.Qbrloss];
       Nest(i).CostSLos=[DistLoadFlowDGSolution.SLosskVA];
end

%% Geting Current best Solution 1
fitness=10^10*ones(User.NumNest,1);
% function
% [BestNest.Cost,Bestnest,Nest.Position,fitness]=get_best_nest(Nest.Position,Nest.Position,fitness);
for j=1:User.NumNest
    
    Nest(j).PLosVolt=feval(User.Function,...
            Nest(j).CostPLos,Nest(j).CostVolt);
        
%     Nest(j).PLosVolt=feval(User.Function,...
%             Nest(j).CostPLos,Nest(j).CostMinVolt,Nest(j).CostQtLos);
    if Nest(j).PLosVolt<=fitness(j)
        fitness(j)=Nest(j).PLosVolt;
        Nest(j).DGnLoc=Nest(j).DGnLoc;
    end
end

% Find the current best
[BestNest.PLosVolt,K]=min(fitness) ;
BestNest.DGnLoc=Nest(K).DGnLoc;
        BestNest.CostPLos       =    Nest(K).CostPLos;
        BestNest.CostPbrLos     =    Nest(K).CostPbrLos;
        BestNest.CostVact       =    Nest(K).CostVact;
        BestNest.CostVolt       =    Nest(K).CostVolt;
        BestNest.CostVSI        =    Nest(K).CostVSI;
        BestNest.CostMinVolt    =    Nest(K).CostMinVolt;
        BestNest.CostVangle     =    Nest(K).CostVangle;
        BestNest.CostQtLos      =    Nest(K).CostQtLos;
        BestNest.CostQbrLos     =    Nest(K).CostQbrLos;
        BestNest.CostSLos       =    Nest(K).CostSLos;
%end of get_best_nest



%% Starting the iteration
for iter =1:User.MaxIter

%% Generate new solutions (but keep the current best)

%fuction New_Nest.Postion=get_cuckoos(Nest.Position,Bestnest.Position,User.Lb,User.Ub);
% Note arg1 is calling in the number os nest (i.e Nest.Postion=User.NumNest)
sigma=(gamma(1+User.beta)*sin(pi*User.beta/2)/(gamma((1+User.beta)/2)*User.beta*2^((User.beta-1)/2)))^(1/User.beta);
for j=1:User.NumNest
    % implementing Levy Flight for each nest
    s=Nest(j).DGnLoc;
    
     % Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/User.beta); 
    
    stepsize=0.01*step.*(s-BestNest.DGnLoc);
    
     s=s+stepsize.*randn(size(s));
     
     %% Application of simple constraints 1
     % function s=simplebounds(s,Lb,Ub)
     % Apply the lower bound
     ns_tmp=s;
     I=ns_tmp<User.Lb;
     ns_tmp(I)=User.Lb(I);
     
     % Apply the upper bounds
     J=ns_tmp>User.Ub;
     ns_tmp(J)=User.Ub(J);
     % Update th is new move
     s=ns_tmp;
     % end of simplebounds
  
  New_Nest(j).DGnLoc=round(s);       % Calling simplebounds
  
  New_DistLoadFlowDGSolution=powerflowDG(New_Nest(j).DGnLoc(1,1),New_Nest(j).DGnLoc(1,2));
    
       New_Nest(j).CostPLos=[New_DistLoadFlowDGSolution.PtLosskW];
       New_Nest(j).CostPbrLos=[New_DistLoadFlowDGSolution.Pbrloss];
       New_Nest(j).CostVact=[New_DistLoadFlowDGSolution.Vactual];
       New_Nest(j).CostVolt=[New_DistLoadFlowDGSolution.VmagPU];
       New_Nest(j).CostVSI=[New_DistLoadFlowDGSolution.VSI];
       New_Nest(j).CostMinVolt=[New_DistLoadFlowDGSolution.minVSI];
       New_Nest(j).CostVangle=[New_DistLoadFlowDGSolution.Vangle];
       New_Nest(j).CostQtLos=[New_DistLoadFlowDGSolution.QtLosskVAr];
       New_Nest(j).CostQbrLos=[New_DistLoadFlowDGSolution.Qbrloss];
       New_Nest(j).CostSLos=[New_DistLoadFlowDGSolution.SLosskVA];
end

%% end of get_cuckoos

%% Geting Current best Solution 2
% calling get_best_nest again... but using New_Nest.Position as input 2nd arguement
for j=1:User.NumNest
%     Nest(j).Cost=feval(User.Function,New_Nest(j).Position);
    
 Nest(j).PLosVolt=feval(User.Function,...
            New_Nest(j).CostPLos,New_Nest(j).CostVolt);
%      Nest(j).PLosVolt=feval(User.Function,...
%             New_Nest(j).CostPLos,New_Nest(j).CostMinVolt,New_Nest(j).CostQtLos);
    
    if Nest(j).PLosVolt<=fitness(j)
        fitness(j)=Nest(j).PLosVolt;
        Nest(j).DGnLoc=New_Nest(j).DGnLoc;
    end
end

% Find the current best
[~,K]=min(fitness) ;
BestNest.DGnLoc=Nest(K).DGnLoc;
BestNest.CostPLos       =    Nest(K).CostPLos;
        BestNest.CostPbrLos     =    Nest(K).CostPbrLos;
        BestNest.CostVact       =    Nest(K).CostVact;
        BestNest.CostVolt       =    Nest(K).CostVolt;
        BestNest.CostVSI        =    Nest(K).CostVSI;
        BestNest.CostMinVolt    =    Nest(K).CostMinVolt;
        BestNest.CostVangle     =    Nest(K).CostVangle;
        BestNest.CostQtLos      =    Nest(K).CostQtLos;
        BestNest.CostQbrLos     =    Nest(K).CostQbrLos;
        BestNest.CostSLos       =    Nest(K).CostSLos;
% End of calling get_best_nest again... but using new_nest as input arguement

%% function new_nest=empty_nests(nest,Lb,Ub,pa) ;
%.............. Discovery and randomization...............

%............ A fraction of worse nests are discovered with a probability pa..............

% Discovered or not -- a status vector
% Converting structure field position to array matrix of NumNest by
% NumbVar
arrayNestDGnLoc=cell2mat({Nest.DGnLoc}.');
K=rand(size(arrayNestDGnLoc))>pa;

% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions.  Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes. 

%% New solution by biased/selective random walks
NumbaofNest=size(arrayNestDGnLoc,1);          %NumbaofNest = User.NumNest
nestn1=arrayNestDGnLoc(randperm(NumbaofNest),:);
nestn2=arrayNestDGnLoc(randperm(NumbaofNest),:);
Nstepsize=rand*(nestn1-nestn2);
new_arrayNestDGnLoc=arrayNestDGnLoc+Nstepsize.*K;

for j=1:size(new_arrayNestDGnLoc,1)
     Ns=new_arrayNestDGnLoc(j,:);
     %Application of simple constraints 2
     % Apply the lower bound 
     Nns_tmp=Ns;
     nI=Nns_tmp<User.Lb;
     Nns_tmp(nI)=User.Lb(nI);
     
     % Apply the upper bounds
     nJ=Nns_tmp>User.Ub;
     Nns_tmp(nJ)=User.Ub(nJ);
     % Update this new move
     Ns=Nns_tmp;
     % end of simplebounds
     New_Nest(j).DGnLoc=round(Ns); 
     
     New_DistLoadFlowDGSolution=powerflowDG(New_Nest(j).DGnLoc(1,1),New_Nest(j).DGnLoc(1,2));
    
       New_Nest(j).CostPLos=[New_DistLoadFlowDGSolution.PtLosskW];
       New_Nest(j).CostPbrLos=[New_DistLoadFlowDGSolution.Pbrloss];
       New_Nest(j).CostVact=[New_DistLoadFlowDGSolution.Vactual];
       New_Nest(j).CostVolt=[New_DistLoadFlowDGSolution.VmagPU];
       New_Nest(j).CostVSI=[New_DistLoadFlowDGSolution.VSI];
       New_Nest(j).CostMinVolt=[New_DistLoadFlowDGSolution.minVSI];
       New_Nest(j).CostVangle=[New_DistLoadFlowDGSolution.Vangle];
       New_Nest(j).CostQtLos=[New_DistLoadFlowDGSolution.QtLosskVAr];
       New_Nest(j).CostQbrLos=[New_DistLoadFlowDGSolution.Qbrloss];
       New_Nest(j).CostSLos=[New_DistLoadFlowDGSolution.SLosskVA];
end

%% Geting Current best Solution 3
% calling get_best_nest again... but using New_Nest.Position as input 2nd arguement
for j=1:User.NumNest
%     New_Nest(j).Cost=feval(User.Function,New_Nest(j).Position);
    
 New_Nest(j).PLosVolt=feval(User.Function,...
            New_Nest(j).CostPLos,New_Nest(j).CostVolt);

%     New_Nest(j).PLosVolt=feval(User.Function,...
%             New_Nest(j).CostPLos,New_Nest(j).CostMinVolt,New_Nest(j).CostQtLos);
    
    if New_Nest(j).PLosVolt<=fitness(j)
        fitness(j)=New_Nest(j).PLosVolt;
        Nest(j).DGnLoc=New_Nest(j).DGnLoc;
    end
end

% Find the current best
[New_BestNest.PLosVolt,K]=min(fitness) ;
New_BestNest.DGnLoc=Nest(K).DGnLoc;

        New_BestNest.CostPLos       =    Nest(K).CostPLos;
        New_BestNest.CostPbrLos     =    Nest(K).CostPbrLos;
        New_BestNest.CostVact       =    Nest(K).CostVact;
        New_BestNest.CostVolt       =    Nest(K).CostVolt;
        New_BestNest.CostVSI        =    Nest(K).CostVSI;
        New_BestNest.CostMinVolt    =    Nest(K).CostMinVolt;
        New_BestNest.CostVangle     =    Nest(K).CostVangle;
        New_BestNest.CostQtLos      =    Nest(K).CostQtLos;
        New_BestNest.CostQbrLos     =    Nest(K).CostQbrLos;
        New_BestNest.CostSLos       =    Nest(K).CostSLos;
% End of calling get_best_nest again... but using new_nest as input arguement
    if New_BestNest.PLosVolt<BestNest.PLosVolt
%         BestNest.PLosVolt=New_BestNest.PLosVolt;
%         BestNest.DGnLoc =New_BestNest.DGnLoc;
        BestNest = New_BestNest;
    end
    
    
    BestCost(iter)=BestNest.PLosVolt;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(iter) ': Best Cost = ' num2str(BestCost(iter))]);
    
end

%Result Display
figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


%% iteration ends
