% visulaization of the k bells trajectories 
% Known BUGS/ Possible Improvements
%
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2022-02-12
% BEGIN USER VAR-------------------------------------------------

addpath('../Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path('../..')
hebec_constants %call the constants function that makes some globals


%%
natoms=100;
v_series=[];
v_series.times=[0,1,2];
inital_v=randn(natoms,2);
inital_v=normr(inital_v);
% check
%sqrt(sum(inital_v.^2,2))
inital_v=inital_v+repmat([0,1],[natoms,1]);
inital_v=cat(1,inital_v,inital_v-repmat([0,2],[natoms,1]));
v_series.vels={inital_v};
up_mask=inital_v(:,2)>0;
v_series.labels=up_mask;



v_second=inital_v;
v_second(up_mask,:)=v_second(up_mask,:)-repmat([0,2],[sum(up_mask),1]);
v_second(~up_mask,:)=v_second(~up_mask,:)+repmat([0,2],[sum(~up_mask),1]);
v_series.vels{2}=v_second;

translate_mask=randn(natoms,2)>0;
v_third=v_second;
v_mix_subset=v_third(translate_mask,:);
up_mask=v_mix_subset(:,2)>0;
yshift=0.01;
v_mix_subset(up_mask,:)=v_mix_subset(up_mask,:)-repmat([yshift,2],[sum(up_mask),1]);
v_mix_subset(~up_mask,:)=v_mix_subset(~up_mask,:)+repmat([yshift,2],[sum(~up_mask),1]);
v_third(translate_mask,:)=v_mix_subset;
v_series.vels{3}=v_third;


%%
t_diff=diff(v_series.times);
elm_distance=cell(1,size(v_series.vels,2)-1);
cum_pos=cell(1,size(v_series.vels,2));
cum_pos{1}=v_series.vels{1}*0;
for ii=1:numel(elm_distance)
    elm_distance{ii}=v_series.vels{ii}*(v_series.times(ii+1)-v_series.times(ii));
    cum_pos{ii+1}=cum_pos{ii}+elm_distance{ii};
end
v_series.cum_pos=cum_pos;
v_series.elm_distance=elm_distance;

%%
stfig('test')
qtimes=linspace(0.1,2.5,200)';
%qtimes=[0.5,1.5,1.9,2.1,2.5];
iimax=numel(qtimes);
for ii=1:iimax
    qtime=qtimes(ii);
    clf
    subplot(1,2,1)
    v_plot=vel_time(v_series,qtime);
    labels_plot=v_series.labels;  %(translate_mask);
    hold on
    plot(v_plot(labels_plot,1),v_plot(labels_plot,2),'ob')
    plot(v_plot(~labels_plot,1),v_plot(~labels_plot,2),'sr')
    axis(gca,'equal')
    hold off
    xlabel('$v_x$')
    ylabel('$v_y$')
    box on
    title('Velocity')
    subplot(1,2,2)
    pos_plot=pos_time(v_series,qtime);
    hold on
    plot(pos_plot(labels_plot,1),pos_plot(labels_plot,2),'ob')
    plot(pos_plot(~labels_plot,1),pos_plot(~labels_plot,2),'sr')
    hold off
    axis(gca,'equal')
    xlabel('$x$')
    ylabel('$y$')
    box on
    title('Position')
    sgtitle(sprintf('$t=%.2f$',qtime))
    pause(0.05)
end


%%




%%

function vel_out=vel_time(v_series,qtime)
idx=find(qtime>v_series.times,1,'last');
vel_out=v_series.vels{idx};
end


function pos_out=pos_time(v_series,qtime)
idx=find(qtime>v_series.times,1,'last');
time_q_delt=qtime-v_series.times(idx);
pos_out=v_series.cum_pos{idx}+vel_time(v_series,qtime)*time_q_delt;
end

