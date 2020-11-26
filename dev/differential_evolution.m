function BestSolHistory = differential_evolution(CostFunction,nVar,int_Position,lb,ub)
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
%% Problem Definition
VarSize=[1 nVar];   % Decision Variables Matrix Size
%% Problem constrains
% int_Position = [0.5000 2.4141 1.0080 -0.8694 15.9905 0.6040 2.1811];%[17.5 1 1.9018 0.6277 0.9000 15 3.6690];%[1 0.9085 0.2019 0.3000 1.4652 0.4849 0.4457];%[4, 0.6,1.06, 0, 9, 0.6,1.06];%[1.8239 0.0595 1.1221 0.1682]; %intial position
% lb = [0.4,0.008,0.03, -0.9, 0.4,0.008,0.03];%[0.1, 0.01, 0.2, -0.9]; %lowerbounds
% ub = [35, 4, 4.0, 0.9, 35, 4, 4.0];%[5, 6, 5, 0.9]; %upperbounds

%% DE Parameters
MaxIt=20;      % Maximum Number of Iterations
nPop=25;        % Population Size
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
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ', Solution = ' num2str(BestSol.Position)]);
    
end
%% Show Results
figure;
%plot(BestCost);
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
