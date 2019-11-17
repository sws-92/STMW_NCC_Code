%% PLOTS for paper
addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))

lnes=lines;

load('20191010_extended_dataset');

% question=input('Would you like to add the 2018-2019 data? Y/N:','s');
% if strcmp(question,'Y')
%     add=load('2019_add_MWdataset');
%     pv_cont_final=[pv_cont_final add.pv_cont_final];
%     low_pv_final=[low_pv_final;add.low_pv_final];
%     heat_content_final=[heat_content_final;add.heat_content_final];
%     stmw_heat_content_final=[stmw_heat_content_final;...
%         add.stmw_heat_content_final];
%     stevens_thickness_final=[stevens_thickness_final;...
%         add.stevens_thickness_final];
%    pv_data_final=[pv_data_final;add.pv_data_final];
%    worthington_intensity_final=[worthington_intensity_final;...
%        add.worthington_intensity_final];
% end

pv_cont_final=pv_cont(:,(time>0));
low_pv_final(low_pv_final(:,1)<0 | low_pv_final(:,1)>2e-10,1)=NaN;
low_pv_final(:,1)=inpaint_nans(low_pv_final(:,1));
worthington_final=worthington(time>0);
heat_content_final(heat_content_final>2.5e7 | heat_content_final<9e6)=NaN;
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
temp_thick_annual=[];
temp_thick_err=[];
pv_thick_annual=[];
pv_thick_err=[];
pv_min_annual=[];
pv_min_err=[];
pv_props_annual=[];
pv_props_err=[];
intensity_annual=[];
intensity_err=[];
heat_content_annual=[];
heat_content_err=[];
stmw_heat_content_annual=[];
stmw_heat_content_err=[];
aou_annual=[];
spice_annual=[];
spice_err=[];

% if strcmp(question,'Y')
     years=1955:2018;
% else
%     years=1955:2018;
% end

for i=years(1):years(end)
        
    % isotherm thickness
    tmp_data=stevens_thickness_final(time_final>=i & time_final<i+1);
    annuali=nanmean(tmp_data);
    temp_thick_annual=[temp_thick_annual;annuali i];
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    temp_thick_err=[temp_thick_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
    % Intensity
    tmp_data=worthington_intensity_final(time_final>i & time_final<i+1);
    intensi=nanmean(tmp_data,1);
    intensity_annual=[intensity_annual;intensi i];
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    intensity_err=[intensity_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals

    % PV Min. 
    tmp_data=low_pv_final((time_final>i & time_final<i+1),1);
    var_annual=nanmean(tmp_data,1);
    pv_min_annual=[pv_min_annual;var_annual i];   
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    pv_min_err=[pv_min_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
    % PV min. properties
    tmp_data=pv_data_final(time_final>i & time_final<i+1,1);
    var_annual=nanmean(tmp_data,1);
    pv_props_annual=[pv_props_annual;var_annual i];   
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    pv_props_err=[pv_props_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
        
    % Heat content
    tmp_data=heat_content_final((time_final>i & time_final<i+1),:);
    var_annual=nanmean(tmp_data);
    heat_content_annual=[heat_content_annual;var_annual i]; 
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    heat_content_err=[heat_content_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
    % STMW heat content 
    tmp_data=stmw_heat_content_final((time_final>i & time_final<i+1),:);
    var_annual=nanmean(stmw_heat_content_final((time_final>i & time_final<i+1),:));
    stmw_heat_content_annual=[stmw_heat_content_annual;var_annual i]; 
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    stmw_heat_content_err=[stmw_heat_content_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
    % Spice 
    tmp_data=s_spice_final((time_final>i & time_final<i+1),:);
    var_annual=nanmean(s_spice_final((time_final>i & time_final<i+1),:));
    spice_annual=[spice_annual;var_annual i]; 
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    spice_err=[spice_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
%     % AOU 
%     var_annual=nanmean(aou_final((time_final>i & time_final<i+1),:));
%     aou_annual=[aou_annual;var_annual i]; 
end

% Load and trim NAO data
load NAO_annual_DJFM.txt
NAO_annual_DJFM(NAO_annual_DJFM(:,1)<1955 | NAO_annual_DJFM(:,1)>2018,:)=[];

%% Weak event- full resolution

% if strcmp(question,'Y')
tmp=[1955:2018]';
weak_mwt=time_final>=1971 & time_final<=1976;
weaker_mwt=time_final>=2014 & time_final<=2018;
% else
% tmp=[1955:2018]';
% weak_mwt=time_final>=1971 & time_final<=1976;
% weaker_mwt=time_final>=2013 & time_final<=2019;
% end
% 

SEM=nanstd(stevens_thickness_final(weak_mwt,1))/sqrt(length(stevens_thickness_final(weak_mwt,1))); % Standard Error
ts=tinv([0.025  0.975],length(stevens_thickness_final(weak_mwt,1))-1); % T-Score
quick_err(1,:)=nanmean(stevens_thickness_final(weak_mwt,1))+ts*SEM; % Confidence Intervals
disp(['Weak 70s event had ',num2str(quick_err(1,1)),'-',num2str(quick_err(1,2)),...
    ' m mean thickness']);

SEM=nanstd(stevens_thickness_final(weaker_mwt,1))/sqrt(length(stevens_thickness_final(weaker_mwt,1))); % Standard Error
ts=tinv([0.025  0.975],length(stevens_thickness_final(weaker_mwt,1))-1); % T-Score
quick_err(2,:)=nanmean(stevens_thickness_final(weaker_mwt,1))+ts*SEM; % Confidence Intervals
disp(['Weak recent event had ',num2str(quick_err(2,1)),'-',num2str(quick_err(2,2)),...
    ' m mean thickness']);

SEM=nanstd(worthington_intensity_final(weak_mwt,1))/sqrt(length(worthington_intensity_final(weak_mwt,1))); % Standard Error
ts=tinv([0.025  0.975],length(worthington_intensity_final(weak_mwt,1))-1); % T-Score
quick_err(3,:)=nanmean(worthington_intensity_final(weak_mwt,1))+ts*SEM; % Confidence Intervals
disp(['Weak 70s event had ',num2str(quick_err(3,1)),'-',num2str(quick_err(3,2)),...
    ' s-1 mean intensity']);

SEM=nanstd(worthington_intensity_final(weaker_mwt,1))/sqrt(length(worthington_intensity_final(weaker_mwt,1))); % Standard Error
ts=tinv([0.025  0.975],length(worthington_intensity_final(weaker_mwt,1))-1); % T-Score
quick_err(4,:)=nanmean(worthington_intensity_final(weaker_mwt,1))+ts*SEM; % Confidence Intervals
disp(['Weak recent event had ',num2str(quick_err(4,1)),'-',num2str(quick_err(4,2)),...
    ' s-1 mean intensity']);

SEM=nanstd(low_pv_final(weak_mwt,1))/sqrt(length(low_pv_final(weak_mwt,1))); % Standard Error
ts=tinv([0.025  0.975],length(low_pv_final(weak_mwt,1))-1); % T-Score
quick_err(5,:)=nanmean(low_pv_final(weak_mwt,1))+ts*SEM; % Confidence Intervals
disp(['Weak 70s event had ',num2str(quick_err(5,1)),'-',num2str(quick_err(5,2)),...
    ' m-1 s-1 mean PV']);

SEM=nanstd(low_pv_final(weaker_mwt,1))/sqrt(length(low_pv_final(weaker_mwt,1))); % Standard Error
ts=tinv([0.025  0.975],length(low_pv_final(weaker_mwt,1))-1); % T-Score
quick_err(6,:)=nanmean(low_pv_final(weaker_mwt,1))+ts*SEM; % Confidence Intervals
disp(['Weak recent event had ',num2str(quick_err(6,1)),'-',num2str(quick_err(6,2)),...
    ' m-1 s-1 mean PV']);

% Load and trim NAO data
load NAO_annual_DJFM.txt
NAO_annual_DJFM(NAO_annual_DJFM(:,1)<1955 | NAO_annual_DJFM(:,1)>2018,:)=[];

%% Plot historical property change
% if strcmp(question,'Y')
    weak_mwt=tmp>=1971 & tmp<=1975;
    weaker_mwt=tmp>=2014 & tmp<=2018;
    recent_mwt=time_final>2010 & time_final<2019;
    tmp_time=1955:2018;
% else
%     weak_mwt=tmp>=1971 & tmp<=1975;
%     weaker_mwt=tmp>=2013 & tmp<=2018;
%     recent_mwt=time_final>1988;
%     tmp_time=1955:2018;
% end
monthi=(1/6);
OCyr=time_final(temp_cont_final(1,:)<19);
% sat_vol=load('20181222_formation_dataset.mat','OC_volume');
% sat_vol=sat_vol.OC_volume-mean(sat_vol.OC_volume);
% sat_vol=interp1(1993:2018,sat_vol,time_final);
% % convert OC volume to Svy if not done already
% if sat_vol(1)>1000
%     sat_vol=sat_vol/(1e6*60*60*24*365);
% end

f1=figure('units','centimeters','outerposition',[0.01 0.01 18 16],'color','w');
axes('position',[0.07 0.94 0.4 0.03])
xtick([]);
ytick([]);
hold on
% contourf([time_final time_final],repmat([0 1],length(time_final),1),...
%     [sat_vol sat_vol],'linestyle','none');
% colormap(m_colmap('diverging'));
for i=1:length(OCyr)
    line([OCyr(i) OCyr(i)],repmat([0 1],length(OCyr),1),'color','k','linewidth',1)
end
box on; xlim([time_final(1) time_final(end)])
text(0.95,0.55,'a)','Units','normalized','FontSize',8,'FontWeight','bold');
axes('position',[0.55 0.94 0.4 0.03])
xtick([]);
ytick([]);
hold on
% contourf([time_final time_final],repmat([0 1],length(time_final),1),...
%     [sat_vol sat_vol],'linestyle','none');
% colormap(m_colmap('diverging'));
for i=1:length(OCyr)
    line([OCyr(i) OCyr(i)],repmat([0 1],length(OCyr),1),'color','k','linewidth',1)
end
box on; xlim([time_final(1) time_final(end)])
text(0.95,0.55,'b)','Units','normalized','FontSize',8,'FontWeight','bold');


axes('position',[0.55 0.74 0.4 0.185]);
hold on
tmp_err=~isnan(pv_props_err(:,2));
jbfill(tmp_time(tmp_err),pv_props_err(tmp_err,2),flipud(pv_props_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),pv_props_err(tmp_err,2),flipud(pv_props_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978,min(pv_props_annual(:,1))+0.02,2,max(pv_props_annual(:,1))-...
    min(pv_props_annual(:,1))],'edgecolor','w','facecolor','w');
hold on

% coeffs=polyfit(time_final(recent_mwt),inpaint_nans(pv_data_final(recent_mwt,1)),1);
[coeffs(1),coeffs(2)]=tsreg(time_final(recent_mwt),...
    inpaint_nans(pv_data_final(recent_mwt,1)),sum(recent_mwt));
fittedy=polyval(coeffs,2010:2018);
l1=plot(2010:2018,fittedy,'-','color',rgb('white'),'linewidth',3); l1=plot(2010:2018,fittedy,'-','color',rgb('pinky red'),'linewidth',2);    
l2=plot(1955:2018,pv_props_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955:2018,pv_props_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--','ytick',16:21);
text(0.95,0.15,'d)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('\Theta (^oC)','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,pv_data_final(:,1)...
    ,2,[0.75 0.75 0.75],'filled');
ylim([15.8 max(pv_data_final(:,1))])
yl=ylim;
yl(1)=yl(1)+(max(yl)/1000);  yl(2)=yl(2)-(max(yl)/1000);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
% for i=1:length(OCyr)
%     f3(i)=patch([OCyr(i)-monthi OCyr(i)+monthi OCyr(i)+monthi OCyr(i)-monthi],...
%         [min(pv_props_err(:,1)) min(pv_props_err(:,1)) max(pv_props_err(:,2)) max(pv_props_err(:,2))],...
%         [0.75 0.75 0.75],'linestyle','none');
% end
uistack(s,'bottom');
% uistack(f3,'bottom')
uistack(f1,'bottom')
uistack(f2,'bottom')
% %uistack(l1,'top')


axes('position',[0.07 0.49 0.4 0.185]);
hold on
tmp_err=~isnan(pv_min_err(:,2));
% jbfill(tmp_time,pv_props_err(tmp_err,2),flipud(pv_props_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),pv_min_err(tmp_err,2),flipud(pv_min_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),pv_min_err(tmp_err,2),flipud(pv_min_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978,min(pv_min_annual(:,1))+2e-12,2,max(pv_min_annual(:,1))-...
    min(pv_min_annual(:,1))],'edgecolor','w','facecolor','w');
hold on

[coeffs(1),coeffs(2)]=tsreg(time_final(recent_mwt),...
    inpaint_nans(low_pv_final(recent_mwt,1)),sum(recent_mwt));
fittedy=polyval(coeffs,2010:2018);
l1=plot(2010:2018,fittedy,'-','color',rgb('white'),'linewidth',3); l1=plot(2010:2018,fittedy,'-','color',rgb('pinky red'),'linewidth',2);    
l2=plot(1955:2018,pv_min_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955:2018,pv_min_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.15,'e)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('Core PV (s^{-1} m^{-1})','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final(low_pv_final(:,1)>0),low_pv_final(low_pv_final(:,1)>0,1)...
    ,2,[0.75 0.75 0.75],'filled');
ylim([-2e-12 max(low_pv_final(:,1))]);
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
% for i=1:length(OCyr)
%     f3(i)=patch([OCyr(i)-monthi OCyr(i)+monthi OCyr(i)+monthi OCyr(i)-monthi],...
%         [min(pv_props_err(:,1)) min(pv_props_err(:,1)) max(pv_props_err(:,2)) max(pv_props_err(:,2))],...
%         [0.75 0.75 0.75],'linestyle','none');
% end
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
%uistack(l1,'top')

axes('position',[0.07 0.24 0.4 0.185]);
hold on
tmp_err=~isnan(intensity_err(:,2));
jbfill(tmp_time(tmp_err),intensity_err(tmp_err,2),flipud(intensity_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),intensity_err(tmp_err,2),flipud(intensity_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978,min(intensity_err(:,1))+2e-9,2,max(intensity_err(:,1))-...
    min(intensity_err(:,1))+2e-9],'edgecolor','w','facecolor','w');
hold on
% [r,p]=Mann_Kendall(intensity_annual(~isnan(intensity_annual(:,1)),1),0.05);
% if r==1
%     coeffs=polyfit(tmp(~isnan(intensity_annual(:,1))),intensity_annual((~isnan(intensity_annual(:,1))),1),1);
%     fittedy=polyval(coeffs,1955:2018);
%     plot([1955:2018]',fittedy,'k--','linewidth',1.5);
% end
% [r,p]=Mann_Kendall(intensity_annual(recent_mwt,1),0.05);

[coeffs(1),coeffs(2)]=tsreg(time_final(recent_mwt),...
    inpaint_nans(worthington_intensity_final(recent_mwt,1)),sum(recent_mwt));
fittedy=polyval(coeffs,2010:2018);
l1=plot(2010:2018,fittedy,'-','color',rgb('white'),'linewidth',3); l1=plot(2010:2018,fittedy,'-','color',rgb('pinky red'),'linewidth',2);    

l2=plot(1955:2018,intensity_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955:2018,intensity_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.9,'g)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('Intensity (s^{-1})','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,worthington_intensity_final(:,1)...
    ,2,[0.75 0.75 0.75],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
% for i=1:length(OCyr)
%     f3(i)=patch([OCyr(i)-monthi OCyr(i)+monthi OCyr(i)+monthi OCyr(i)-monthi],...
%         [min(pv_props_err(:,1)) min(pv_props_err(:,1)) max(pv_props_err(:,2)) max(pv_props_err(:,2))],...
%         [0.75 0.75 0.75],'linestyle','none');
% end
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
%uistack(l1,'top')

axes('position',[0.07 0.74 0.4 0.185]);
hold on
tmp_err=~isnan(temp_thick_err(:,2));
jbfill(tmp_time(tmp_err),temp_thick_err(tmp_err,2),flipud(temp_thick_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),temp_thick_err(tmp_err,2),flipud(temp_thick_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978,min(temp_thick_err(:,1))+50,2,max(temp_thick_err(:,1))-...
    50],'edgecolor','w','facecolor','w');
hold on

[coeffs(1),coeffs(2)]=tsreg(time_final(recent_mwt),...
    inpaint_nans(stevens_thickness_final(recent_mwt,1)),sum(recent_mwt));
fittedy=polyval(coeffs,2010:2018);
l1=plot(2010:2018,fittedy,'-','color',rgb('white'),'linewidth',3); l1=plot(2010:2018,fittedy,'-','color',rgb('pinky red'),'linewidth',2);    

l2=plot(1955:2018,temp_thick_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955:2018,temp_thick_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'gridlinestyle','--','ytick',0:100:500);
text(0.95,0.9,'c)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('Thickness (m)','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,stevens_thickness_final(:,1)...
    ,2,[0.75 0.75 0.75],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
% for i=1:length(OCyr)
%     f3(i)=patch([OCyr(i)-monthi OCyr(i)+monthi OCyr(i)+monthi OCyr(i)-monthi],...
%         [min(pv_props_err(:,1)) min(pv_props_err(:,1)) max(pv_props_err(:,2)) max(pv_props_err(:,2))],...
%         [0.75 0.75 0.75],'linestyle','none');
% end
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
%uistack(l1,'top')


axes('position',[0.55 0.49 0.4 0.185]);
hold on
tmp_err=~isnan(stmw_heat_content_err(:,2));
jbfill(tmp_time(tmp_err),stmw_heat_content_err(tmp_err,2),flipud(stmw_heat_content_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),stmw_heat_content_err(tmp_err,2),flipud(stmw_heat_content_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978,min(stmw_heat_content_err(:,1))+1e6,2,max(stmw_heat_content_err(:,1))-...
    min(stmw_heat_content_err(:,1))],'edgecolor','w','facecolor','w');
hold on
[coeffs(1),coeffs(2)]=tsreg(time_final(recent_mwt),...
    inpaint_nans(stmw_heat_content_final(recent_mwt,1)),sum(recent_mwt));
fittedy=polyval(coeffs,2010:2018);
l1=plot(2010:2018,fittedy,'-','color',rgb('white'),'linewidth',3); l1=plot(2010:2018,fittedy,'-','color',rgb('pinky red'),'linewidth',2);  
l2=plot(1955:2018,stmw_heat_content_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955:2018,stmw_heat_content_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.9,'f)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel({'{\ith_{STMW}} (J kg^{-1} m)'},'FontSize',8,'FontWeight','bold',...
    'interpreter','tex');
set(gca,'FontSize',8);
s=scatter(time_final,stmw_heat_content_final(:,1)...
    ,2,[0.75 0.75 0.75],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
% for i=1:length(OCyr)
%     f3(i)=patch([OCyr(i)-monthi OCyr(i)+monthi OCyr(i)+monthi OCyr(i)-monthi],...
%         [min(pv_props_err(:,1)) min(pv_props_err(:,1)) max(pv_props_err(:,2)) max(pv_props_err(:,2))],...
%         [0.75 0.75 0.75],'linestyle','none');
% end
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
%uistack(l1,'top')

axl=axes('position',[0.55 0.24 0.4 0.185]);
hold on
tmp_err=~isnan(heat_content_err(:,2));
jbfill(tmp_time(tmp_err),heat_content_err(tmp_err,2),flipud(heat_content_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),heat_content_err(tmp_err,2),flipud(heat_content_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978,min(heat_content_err(:,1))+5e4,2,max(heat_content_err(:,1))-...
    min(heat_content_err(:,1))],'edgecolor','w','facecolor','w');
hold on
[coeffs(1),coeffs(2)]=tsreg(time_final(recent_mwt),...
    inpaint_nans(heat_content_final(recent_mwt,1)),sum(recent_mwt));
fittedy=polyval(coeffs,2010:2018);
l1=plot(2010:2018,fittedy,'-','color',rgb('white'),'linewidth',3); l1=plot(2010:2018,fittedy,'-','color',rgb('pinky red'),'linewidth',2);  
l2=plot(1955:2018,heat_content_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955:2018,heat_content_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.15,'h)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel({'{\ith_{s}} (J kg^{-1} m)'},'FontSize',8,'FontWeight','bold',...
    'interpreter','tex');
set(gca,'FontSize',8);
s=scatter(time_final,heat_content_final(:,1)...
    ,2,[0.75 0.75 0.75],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/500);  yl(2)=yl(2)-(max(yl)/500);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
% for i=1:length(OCyr)
%     f3(i)=patch([OCyr(i)-monthi OCyr(i)+monthi OCyr(i)+monthi OCyr(i)-monthi],...
%         [min(heat_content_err(:,1)) min(heat_content_err(:,1)) max(heat_content_err(:,2)) max(heat_content_err(:,2))],... 
%         [0.75 0.75 0.75],'linestyle','none');
% end
uistack(s,'bottom');
% uistack(f3,'bottom')
uistack(f1,'bottom')
uistack(f2,'bottom')
%uistack(l1,'top')

axes('position',[0.3 0.04 0.4 0.14]);
hold on
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [min(NAO_annual_DJFM(:,2))+0.15 min(NAO_annual_DJFM(:,2))+0.15 max(NAO_annual_DJFM(:,2))-0.1 max(NAO_annual_DJFM(:,2))-0.1],...
    rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [min(NAO_annual_DJFM(:,2))+0.15 min(NAO_annual_DJFM(:,2))+0.15  max(NAO_annual_DJFM(:,2))-0.1 max(NAO_annual_DJFM(:,2))-0.1],...
    rgb('very light blue'),'linestyle','none');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.05,0.9,'i)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel({'NAO Index';'(DJFM)'},'FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
% Create high res NAO signal for fill
timei=linspace(NAO_annual_DJFM(1,1)+0.5,NAO_annual_DJFM(end,1)+0.5,1000);
NAOi=interp1(NAO_annual_DJFM(:,1)+0.5,NAO_annual_DJFM(:,2),timei');
posneg=NAOi>0;
up=zeros(length(posneg),1);
low=zeros(length(posneg),1);
up(posneg)=NAOi(posneg);
low(~posneg)=NAOi(~posneg);
[ph,msg]=jbfill(timei,up,zeros(length(posneg),1),'white','white',1,1);
[ph,msg]=jbfill(timei,low,zeros(length(posneg),1),'white','white',1,1);
[ph,msg]=jbfill(timei,up,zeros(length(posneg),1),'red','red',1,0.5);
[ph,msg]=jbfill(timei,low,zeros(length(posneg),1),'blue','blue',1,0.5);
box on
xlim([1955 2018]);

%% 

export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\FigureTS.eps -deps
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\FigureTS.pdf -dpdf

%% Spice
figure;
hold on
tmp_err=~isnan(spice_err(:,2));
jbfill(tmp_time(tmp_err),spice_err(tmp_err,2),flipud(spice_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),spice_err(tmp_err,2),flipud(spice_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978,min(spice_err(:,1)),2,max(spice_err(:,1))-...
    min(spice_err(:,1))],'edgecolor','w','facecolor','w');
hold on
coeffs=polyfit(time_final(recent_mwt),inpaint_nans(s_spice_final(recent_mwt,1)),1);
fittedy=polyval(coeffs,2011:2017);
% l1=plot(2011:2017,fittedy,'-','color',rgb('bright red'),'linewidth',2);
l2=plot(1955:2018,spice_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955:2018,spice_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.15,'h)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel({'Spiciness'},'FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,s_spice_final(:,1)...
    ,2,[0.75 0.75 0.75],'filled');
yl=ylim;
% yl(1)=yl(1)+(max(yl));  yl(2)=yl(2);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
% for i=1:length(OCyr)
%     f3(i)=patch([OCyr(i)-monthi OCyr(i)+monthi OCyr(i)+monthi OCyr(i)-monthi],...
%         [min(spice_err(:,1)) min(spice_err(:,1)) max(spice_err(:,2)) max(spice_err(:,2))],... 
%         [0.75 0.75 0.75],'linestyle','none');
% end
uistack(s,'bottom');
% uistack(f3,'bottom')
uistack(f1,'bottom')
uistack(f2,'bottom')
%uistack(l1,'top')
