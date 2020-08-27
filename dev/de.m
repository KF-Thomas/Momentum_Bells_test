%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA107
% Project Title: Implementation of Differential Evolution (DE) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
% clc;
% clear;
% close all;
%% Optimization Options
cost_opts.diffraction_order = 9;
cost_opts.orders = [-2,-1];%[-2,-1,0]; %the difraction orders we want top optimize transfer to
%% Problem Definition
CostFunction=@(b) Raman_Nath_Cost(b,cost_opts);    % Cost Function
nVar=4;            % Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size
%% Problem constrains
int_Position = [1.8239 0.0595 1.1221 0.1682]; %intial position
lb = [0.1, 0.01, 0.2, -0.9]; %lowerbounds
ub = [5, 6, 5, 0.9]; %upperbounds

%% DE Parameters
MaxIt=1000;      % Maximum Number of Iterations
nPop=55;        % Population Size
beta_min=0.1;   % Lower Bound of Scaling Factor
beta_max=0.9;   % Upper Bound of Scaling Factor
pCR=0.2;        % Crossover Probability
%% Initialization
empty_individual.Position=[];
empty_individual.Cost=[];
BestSol.Cost=inf;
BestSolHistory={};
pop=repmat(empty_individual,nPop,1);
for i=1:nPop
    pop(i).Position=unifrnd(lb,ub,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
        BestSolHistory=[BestSolHistory;BestSol];
    end
    
end
BestCost=zeros(MaxIt,1);
%% DE Main Loop
for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        %beta=unifrnd(beta_min,beta_max);
        beta=unifrnd(beta_min,beta_max,VarSize);
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        % intial point
        if it == 1 && i == 1 && ~isempty(int_Position)
            z = int_Position;
        end
        
        % Check bounds
        for j=1:numel(z)
            if z(j)<lb(j)
                z(j) = lb(j);
            elseif z(j)>ub(j)
                z(j) = ub(j);
            end
        end
        
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
               BestSolHistory=[BestSolHistory;BestSol];
            end
        end
        
    end
    
    % Update Best Cost
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end
%% Show Results
figure;
%plot(BestCost);
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
