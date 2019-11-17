%% Different STMW definitions
addpath(genpath('C:\Users\samst\Documents'))
clear

load 20191010_extended_dataset

pv_cont_final=pv_cont(:,(time>0));
low_pv_final(low_pv_final(:,1)<0 |...
    low_pv_final(:,1)>2e-10,1)=NaN;
low_pv_final(:,1)=inpaint_nans(low_pv_final(:,1));
worthington_final=worthington(time>0);
heat_content_final(heat_content_final>2.5e7 |...
    heat_content_final<9e6)=NaN;
heat_content_final=inpaint_nans(heat_content_final);
stmw_heat_content_final(heat_content_final<9e6)=NaN;
stmw_heat_content_final=inpaint_nans(stmw_heat_content_final);

pv_data_final((pv_data_final(:,1)==0),1)=NaN;
pv_data_final(:,1)=inpaint_nans(pv_data_final(:,1));
pv_data_final((pv_data_final(:,2)==0),2)=NaN;
pv_data_final(:,2)=inpaint_nans(pv_data_final(:,2));
pv_data_final((pv_data_final(:,3)==0),3)=NaN;
pv_data_final(:,3)=inpaint_nans(pv_data_final(:,3));
pv_data_final((pv_data_final(:,4)==0),4)=NaN;
pv_data_final(:,4)=inpaint_nans(pv_data_final(:,4));
pv_data_final((pv_data_final(:,5)==0),5)=NaN;
pv_data_final(:,5)=inpaint_nans(pv_data_final(:,5));

%% Annual datasets
% if strcmp(question,'Y')
years=1955:2018;
% else
%     years=1955:2018;
% end

temp_thick_annual=[];
temp_thick_err=[];
pv_min_annual=[];
pv_min_err=[];
t_def_annual=[];
pv_props_annual=[];
pv_props_err=[];
intensity_annual=[];
intensity_err=[];

for i=years(1):years(end)
    
    % isotherm thickness
    tmp_data=thickness_def_final(time_final>=i &  time_final<i+1,:);
    tmp_data=[tmp_data stevens_thickness_final(time_final>=i &  time_final<i+1)];
    annuali=mean(tmp_data,1,'omitnan');
    temp_thick_annual=[temp_thick_annual;annuali];

    % PV Min.
    tmp_data=pv_def_final(( time_final>i &  time_final<i+1),:);
    tmp_data=[tmp_data low_pv_final(time_final>=i &  time_final<i+1,1)];
    var_annual=mean(tmp_data,1,'omitnan');
    pv_min_annual=[pv_min_annual;var_annual];
    
     % Core temp.
    tmp_data=t_def_final(( time_final>i &  time_final<i+1),:);
    tmp_data=[tmp_data pv_data_final(time_final>=i &  time_final<i+1,1)];
    var_annual=mean(tmp_data,1,'omitnan');
    t_def_annual=[t_def_annual;var_annual];
end

%% tmp=1955:2018;
tmp=1955:2018;
weak_mwt=tmp>=1971 & tmp<=1975;
weaker_mwt=tmp>=2014 & tmp<=2018;

figure('units','centimeters','outerposition',[0 0 12 12]);
set(gcf,'color','w');
cm=linspecer(size(temp_thick_annual,2));

ax(1)=subplot(3,1,1);
hold on
ax(1).XAxis.Visible='off';
for i=1:size(temp_thick_annual,2)
    l(i)=plot(1955:2018,movmean(...
        temp_thick_annual(:,i),3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(i,:));
end
axis tight
yl=ylim;
% ylim([0 max(yl)]);
ax(1).Color='none';
ylabel('Thickness (m)','FontSize',8,'FontWeight','bold');
ax(1).YGrid='on';
ax(1).GridAlpha = 0.06;
% xtick([0 100 200]);
% legend('location','best');

ax(2)=subplot(3,1,2);
hold on
ax(2).XAxis.Visible='off';
for i=1:size(temp_thick_annual,2)
    l(i)=plot(1955:2018,movmean(...
        pv_min_annual(:,i),3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(i,:));
end
axis tight
ax(2).Color='none';
ylabel({'Core PV';'(s^{-1} m^{-1})'},'FontSize',8,'FontWeight','bold');
ax(2).YGrid='on';
ax(2).GridAlpha = 0.06;

ax(3)=subplot(3,1,3);
hold on
ax(3).XAxis.Visible='off';
for i=1:size(temp_thick_annual,2)
    l(i)=plot(1955:2018,movmean(...
        t_def_annual(:,i),3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(i,:));
end
axis tight
ax(3).Color='none';
ylabel({'\theta (^\circC)'},'FontSize',8,'FontWeight','bold');
ax(3).YGrid='on';
ax(3).GridAlpha = 0.06;

iax=axes;
plot(1955:2018,NaN(length(1955:2018),1));
yl=ylim;
yl(1)=yl(1)+(max(yl)/1000);  yl(2)=yl(2)-(max(yl)/1000);
f1=patch([min(tmp(weak_mwt)) max(tmp(weak_mwt)) max(tmp(weak_mwt)) min(tmp(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp(weaker_mwt)) max(tmp(weaker_mwt)) max(tmp(weaker_mwt)) min(tmp(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
linkaxes([iax ax],'x');
iax.YAxis.Visible='off';
iax.XGrid='on';
iax.XMinorGrid='on';
iax.MinorGridLineStyle='-';
iax.GridAlpha = 0.075;
iax.MinorGridAlpha = 0.04;
% y=legend(l,{'BATS/HS','Argo','EN4g10','ORAS4'},'FontSize',6,'Location','best');
% title(y,'Dataset');
uistack(iax,'bottom');
box off


%% Find means for both periods
weak_mwt=time_final>=1971 & time_final<=1975;
weaker_mwt=time_final>=2014 & time_final<2019;

tmp_thick=[thickness_def_final(weak_mwt,:) stevens_thickness_final(weak_mwt)];
tmp_pv=[pv_def_final(weak_mwt,:) low_pv_final(weak_mwt,1)];

SEM=std(tmp_thick,1,'omitnan')/sqrt(size(tmp_thick,1)); % Standard Error
SEMp=std(tmp_pv,1,'omitnan')/sqrt(size(tmp_pv,1)); % Standard Error
ts=tinv([0.05  0.95],size(tmp_thick,1)-1); % T-Score
for i=1:size(tmp_thick,2)
    thick_err_weak(i,:)=mean(tmp_thick(:,i),1,'omitnan')+ts*SEM(i); % Confidence Intervals
    pv_err_weak(i,:)=mean(tmp_pv(:,i),1,'omitnan')+ts*SEMp(i); % Confidence Intervals
end

tmp_thick=[thickness_def_final(weaker_mwt,:) stevens_thickness_final(weaker_mwt)];
tmp_pv=[pv_def_final(weaker_mwt,:) low_pv_final(weaker_mwt,1)];

SEM=std(tmp_thick,1,'omitnan')/sqrt(size(tmp_thick,1)); % Standard Error
SEMp=std(tmp_pv,1,'omitnan')/sqrt(size(tmp_pv,1)); % Standard Error
ts=tinv([0.05  0.95],size(tmp_thick,1)-1); % T-Score
for i=1:size(tmp_thick,2)
    thick_err_weaker(i,:)=mean(tmp_thick(:,i),1,'omitnan')+ts*SEM(i); % Confidence Intervals
    pv_err_weaker(i,:)=mean(tmp_pv(:,i),1,'omitnan')+ts*SEMp(i); % Confidence Intervals
end

disp('All PV in 2010s more than 1970s (95%)');
disp('3/5 PV in 2010s less than 1970s (95%, methods 3,4,5)');
disp('2/5 PV in 2010s less than 1970s (80%, methods 1,2)');

%% 

for i=years(1)+2:years(end)-2
    
    % isotherm thickness
    tmp_data=thickness_def_final(time_final>=i-2.5 &  time_final<i+2.5,:);
    annuali=mean(tmp_data,1,'omitnan');
    temp_thick_annual=[temp_thick_annual;annuali];

    % PV Min.
    tmp_data=pv_def_final(( time_final>i-2.5 &  time_final<i+2.5),:);
    var_annual=mean(tmp_data,1,'omitnan');
    pv_min_annual=[pv_min_annual;var_annual];
end
