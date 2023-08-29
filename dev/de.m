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
% cost_opts.diffraction_order = 9;
% cost_opts.orders = [-2,-1];%[-2,-1,0]; %the difraction orders we want top optimize transfer to
% cost_opts.goal = 'transfer';%'mirror';%'splitting';
% %% Problem Definition
% cost_opts.diffraction_order = 9;
% cost_opts.orders = [-1,-2];
% cost_opts.goal = 'mirror';%'transfer';%'splitting';

CostFunction=@(b) S_opt(b,out_data);%Raman_Nath_Cost(b,cost_opts);    % Cost Function
nVar=11;            % Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size
%% Problem constrains
int_Position = [-1 600 1e3 9e3 0 0.1 1 1.5 1 1 1];%[17.5 1 1.9018 0.6277 0.9000 15 3.6690];%[1 0.9085 0.2019 0.3000 1.4652 0.4849 0.4457];%[4, 0.6,1.06, 0, 9, 0.6,1.06];%[1.8239 0.0595 1.1221 0.1682]; %intial position
lb = [-1,    1,   1, 1.5e3, 0,   0,  0.1, 1, 0.1, 0.1, 0.1];%[0.1, 0.01, 0.2, -0.9]; %lowerbounds
ub = [300, 1e6, 5e3,   1e6, pi, pi,  3.5, 2.5, 15, 15, 15];%[5, 6, 5, 0.9]; %upperbounds

%% DE Parameters
MaxIt=100;      % Maximum Number of Iterations
nPop=55;        % Population Size
beta_min=0.2;   % Lower Bound of Scaling Factor
beta_max=0.8;   % Upper Bound of Scaling Factor
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
plot(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
%%
S_opt(BestVals(10,:),out_data)

%%
function cost = S_opt(b,out_data)

corr_opts.verbose = false;
    corr_opts.print_update = false;
    corr_opts.timer=false;
    corr_opts.plots = false;
    corr_opts.fit = false;

corr_opts.rad_smoothing=nan;
    corr_opts.one_d_smoothing=nan;

    corr_opts.direction_labels = {'z','x','y'};
    corr_opts.low_mem=true;

    corr_opts.do_pre_mask=false;
    corr_opts.sorted_dir=nan;
    corr_opts.sort_norm=0;
    

    % variables for calculating the error
    corr_opts.calc_err = false;
    corr_opts.samp_frac_lims=[0.65,0.9];
    corr_opts.num_samp_frac=5;
    corr_opts.num_samp_rep=5;
    
    corr_opts.attenuate_counts=1;
    corr_opts.one_d_edges = linspace(-0.015,0.015,15.*2);
    corr_opts.one_d_dimension = 2;
    corr_opts.type='1d_vol_bb';%'1d_cart_bb';%'radial_bb';%

%radial mask

%centering mask

% Halo N lims
halo_N_lims = [b(1),b(2)]; %2

% Num lims
data_N_lims = [b(3),b(4)];

%theta bins
corr_opts.theta_bins = [b(5),b(6)]; %2
corr_opts.theta_mask_polarity = 0;

%phi bins
L_ang = b(7); % in units of 10 degrees
Ang_cen = b(8); % in units of 10 degrees
phi_bins_S = [-L_ang,L_ang;Ang_cen-L_ang,Ang_cen+L_ang;-L_ang,L_ang].*pi/18; %2

%bin size
opts_E.delta_kd =  [b(9)*1e-3,b(10)*1e-3,b(11)*1e-3]; %3

    dkx = opts_E.delta_kd(2);
    dky = opts_E.delta_kd(3);
    dkz = opts_E.delta_kd(1);
    
    
    corr_opts.bin_lims = 15;

    corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];

    corr_opts.sample_proportion=0.0001;%1.0;%0.65;%1500;
    corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')

S_indx = [4,5,6];

clear sph_S_corrs
for ii = 1:length(S_indx)

N_mask = out_data{S_indx(ii)}.bottom_halo.num_counts<halo_N_lims(2) & ...
out_data{S_indx(ii)}.top_halo.num_counts<halo_N_lims(2) & ...
out_data{S_indx(ii)}.bottom_halo.num_counts>halo_N_lims(1) & ...
out_data{S_indx(ii)}.top_halo.num_counts>halo_N_lims(1);

atom_num = cellfun(@(x) size(x,1), out_data{S_indx(ii)}.top_halo.counts_txy);
atom_mask = atom_num>data_N_lims(1) & atom_num<data_N_lims(2);

corr_opts.phi_bins = phi_bins_S(ii,:);

current_data.top_halo = struct_mask(out_data{S_indx(ii)}.top_halo,N_mask&atom_mask);
current_data.bottom_halo = struct_mask(out_data{S_indx(ii)}.bottom_halo,N_mask&atom_mask);

if isempty(current_data.top_halo.counts_txy)
     cost = 5e3;
     return
end

sph_S_corrs{ii}.top = spherical_section_g2(corr_opts,current_data.top_halo.counts_vel');
sph_S_corrs{ii}.btm = spherical_section_g2(corr_opts,current_data.bottom_halo.counts_vel');
both_halo_counts = [current_data.top_halo.counts_vel';current_data.bottom_halo.counts_vel'];
sph_S_corrs{ii}.btw = spherical_section_g2(corr_opts,both_halo_counts);
end

g2_S_top = cell2mat(cellfun(@(x) x.top.data.in_shot_corr.one_d_corr_density(:),sph_S_corrs,'UniformOutput',false));
g2_S_btm = cell2mat(cellfun(@(x) x.btm.data.in_shot_corr.one_d_corr_density(:),sph_S_corrs,'UniformOutput',false));
g2_S_btw = cell2mat(cellfun(@(x) x.btw.data.in_shot_corr.one_d_corr_density(:),sph_S_corrs,'UniformOutput',false));

SG=g2_S_top+g2_S_btm;
XG=2.*g2_S_btw;

N_shots_top = cell2mat(cellfun(@(y) cellfun(@(x) size(x,1),y.top.data.in_shot_corr.shots).',sph_S_corrs,'UniformOutput',false));
N_shots_btm = cell2mat(cellfun(@(y) cellfun(@(x) size(x,1),y.btm.data.in_shot_corr.shots).',sph_S_corrs,'UniformOutput',false));
N_shots_btw = cell2mat(cellfun(@(y) cellfun(@(x) size(x,1),y.btw.data.in_shot_corr.shots).',sph_S_corrs,'UniformOutput',false));

SG_unc = sqrt(g2_S_top./N_shots_top+g2_S_btm./N_shots_btm);
XG_unc = 2.*sqrt(g2_S_btw./N_shots_btw);

E_sph = ((SG-XG)./(SG+XG));
E_unc = E_sph.*sqrt((SG_unc.^2+XG_unc.^2)./(SG-XG).^2+(SG_unc.^2+XG_unc.^2)./(SG+XG).^2);

S_sph = sum(E_sph.*[-1,2,1],2);
S_unc = sqrt(sum((E_unc.*[1,2,1]).^2,2));
cost = min((2-S_sph)./S_unc);
end