%% This script creates plots, statistics, tables etc. 
% SWS, RJJ, NRB. STMW PAPER
% This script is a compilation of all of the code written to produce the
% figures and statistics for Stevens et al. 2018. Email SWS at
% sstevens@eoas.ubc.ca if you have any questions. 

% This code was lifted from several different scripts, and was
% developed episodically over many years, thus it is not well commented.
% The dataset scripts included with this commit create the datasets to
% generate the items below. 

addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))

%% FIGURE - Map
clear

load 20190504_formation_dataset.mat

% convert OC volume to Svy if not done already
if OC_volume(1)>1000
    OC_volume=OC_volume/(1e6*60*60*24*365);
    vol_err=vol_err/(1e6*60*60*24*365);
end

% Set up figure
f1=figure('units','centimeters','outerposition',[0.01 0.01 18 13]);
colormap(cmocean('Thermal'));
lon_lim=[-78 -42];
lat_lim=[25.5 45];
m_proj('lambert','lon',lon_lim,'lat',lat_lim,'rect','on');   % Projection
hold on
RGB=rgb('light grey');

% Find formation box idxs
lon_idx=sat_lon>-75 & sat_lon<-45;
lat_idx=sat_lat>30 & sat_lat<41.5;

% Contour formation zones
[C,h]=m_contourf(sat_lon,sat_lat,mean_sat_sst,'LineColor','w',...
    'levelstep',1,'LevelList',2:1:23,'LineColor','none');
[ax,h]=m_contfbar(-0.05,[0.2 0.8],C,h,'axfrac',.02,'fontsize',8,'fontweight','bold');
ylabel(ax,{'SST';'(^oC)'},'rotation',0,'color','k');
[C,c]=m_contour(sat_lon(lon_idx),sat_lat(lat_idx),mean_sat_sst(lat_idx,lon_idx),[17 19],'w--','linewidth',1.5);
[B,b]=m_contour(sat_lon(lon_idx),sat_lat(lat_idx),mean_sat_sst(lat_idx,lon_idx),[18 18],'w-.','linewidth',1.5);
[F,f]=m_contour(sat_lon(lon_idx),sat_lat(lat_idx),sst_mean_pentad(lat_idx,lon_idx),[17 19],'k--','linewidth',1);
[F,f]=m_contour(sat_lon(lon_idx),sat_lat(lat_idx),sst_mean_pentad(lat_idx,lon_idx),[18 18],'k-.','linewidth',1);

m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','on');
[x,y]=meshgrid(sat_lon,sat_lat);
m_line([x(1,:) x(:,end)' fliplr(x(1,:)) x(:,1)'],...
    [y(end,:) fliplr(y(:,end)') y(1,:) y(:,end)'],...
    'color','k','linestyle','--','linewidth',1.5);
[x,y]=meshgrid(-75:-45,30:41);
m_line([x(1,:) x(:,end)' fliplr(x(1,:)) x(:,1)'],...
    [y(end,:) fliplr(y(:,end)') y(1,:) y(:,end)'],...
    'color','k','linestyle','-.','linewidth',1.5);
m_scatter(-64.16666,31.666,25,'k','filled');
m_text(-64.16666+0.5,31.666-0.1,'BATS','color','k','fontweight','bold','fontsize',8);
m_scatter(-64.5,32.16666,25,'k','filled');
m_text(-64.5+0.5,32.16666+0.1,'Hydrostation ''S''','color','k','fontweight','bold','fontsize',8);
% text(0.9,0.9,'a)','Units','normalized','FontSize',8,'FontWeight','bold');

% Add BV stations
BVal=load('BV_trend_data.mat');
m_scatter(BVal.BVintens_trend(:,3),BVal.BVintens_trend(:,4),15,...
    'k','filled');
m_text(BVal.BVintens_trend([1 3:7],3)+0.3,BVal.BVintens_trend([1 3:7],4),...
    {'BV1','BV3','BV4','BV5','BV6','BV7'},...
    'color','k','fontsize',8);

%% FIGURE- PV profiles
clear

if ~exist('HS','var')
    HS=load('20190504_extended_dataset');
end
timeyear=1955:2017;
[hc,edges]=histcounts(floor(HS.time_final),1955:2017);

f1=figure('units','centimeters','outerposition',[0.01 0.01 9 15]);
subplot(3,1,1)
histogram(floor(HS.time_final),1955:2017,'facecolor',rgb('grey'));
set(gca,'yscale','log','ytick',[1 10 50 100 500]);
hold on
% Have to add these in as matlab doesnt plot them...
patch([1974 1975 1975 1974],[0 0 hc(edges==1974) hc(edges==1974)],...
    rgb('grey'),'edgecolor','k','linewidth',1);
patch([1981 1982 1982 1981],[0 0 hc(edges==1981) hc(edges==1981)],...
    rgb('grey'),'edgecolor','k','linewidth',1);
axis tight
ylim([0 600]);
ylabel('Number of profiles','FontSize',8,'FontWeight','bold');
text(0.02,0.9,'a)','Units','normalized','FontSize',8,'FontWeight','bold')
% xlabel('Year','FontSize',8,'FontWeight','bold');
grid on
set(gca,'GridLineStyle','--','tickdir','out');

mn_profile=nanmean(HS.pv_cont_final,2);
pent_profile=nanmean(HS.pv_cont_final(:,HS.time_final>2013 &...
    HS.time_final<2018),2);
dpth_prf=5:5:600;
idxf=find(mn_profile<1e-10,1,'first');
idxl=find(mn_profile<1e-10,1,'last');
subplot(3,1,2:3)
hold on
xlim([0 5e-10]);
p1=patch([xlim fliplr(xlim)],[dpth_prf(idxf) dpth_prf(idxf)...
    dpth_prf(idxl) dpth_prf(idxl)],rgb('light blue'),'edgecolor','none',...
    'facealpha',0.5');
plot(HS.pv_cont_final,dpth_prf,':','color',[0.5 0.5 0.5],'linewidth',0.5);
plot(mn_profile,dpth_prf,'w','linewidth',3);
%plot(mn_profile(1:idxf),dpth_prf(1:idxf),'r','linewidth',2);
plot(mn_profile(1:idxf),dpth_prf(1:idxf),'b:','linewidth',1.5);
plot(mn_profile(idxf:idxl),dpth_prf(idxf:idxl),'b','linewidth',2);
%plot(mn_profile(idxl:end),dpth_prf(idxl:end),'r','linewidth',2);
plot(mn_profile(idxl:end),dpth_prf(idxl:end),'b:','linewidth',1.5);

idxf=find(pent_profile<1e-10,1,'first');
idxl=find(pent_profile<1e-10,1,'last');
p2=patch([xlim fliplr(xlim)],[dpth_prf(idxf) dpth_prf(idxf)...
    dpth_prf(idxl) dpth_prf(idxl)],'w','edgecolor','none');
p3=patch([xlim fliplr(xlim)],[dpth_prf(idxf) dpth_prf(idxf)...
    dpth_prf(idxl) dpth_prf(idxl)],rgb('light red'),'edgecolor','none',...
    'facealpha',0.5');
plot(pent_profile,dpth_prf,'w','linewidth',3);
%plot(pent_profile(1:idxf),dpth_prf(1:idxf),'color',lnes(5,:),'linewidth',2);
plot(pent_profile(1:idxf),dpth_prf(1:idxf),'r:','linewidth',1.5);
plot(pent_profile(idxf:idxl),dpth_prf(idxf:idxl),'r','linewidth',2);
%plot(pent_profile(idxl:end),dpth_prf(idxl:end),'color',lnes(5,:),'linewidth',2);
plot(pent_profile(idxl:end),dpth_prf(idxl:end),'r:','linewidth',1.5);
set(gca,'ydir','reverse');
axis tight
ylim([100 600]);
xlim([0 5e-10]);
text(0.95,0.05,'b)','Units','normalized','FontSize',8,'FontWeight','bold')
ylabel('Depth (m)','FontSize',8,'FontWeight','bold');
xlabel('PV (m^{-1} s^{-1})','FontSize',8,'FontWeight','bold');
uistack(p1,'bottom');
uistack(p2,'bottom');uistack(p2,'up');
uistack(p3,'bottom');uistack(p3,'up',2);
% text(3.5e-10,275,'1955-2017','color','b','fontweight','bold','fontsize',8);
annotation('textbox',[0.65 0.46 0.01 0.01],'String','1955-2017','Fitboxtotext','on',...
    'color','b','fontweight','bold','fontsize',8,'backgroundcolor','w');
annotation('textbox',[0.65 0.34 0.01 0.01],'String','2013-2017','Fitboxtotext','on',...
    'color','r','fontweight','bold','fontsize',8,'backgroundcolor','w');
box on
set(gca,'tickdir','out','fontsize',8);

%% FIGURE- Time-Series
clear
lnes=lines;

load('20190504_extended_dataset');

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

% Annual datasets
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

for i=1955:2017
        
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
    
    % AOU 
    var_annual=nanmean(aou_final((time_final>i & time_final<i+1),:));
    aou_annual=[aou_annual;var_annual i]; 
end

% Load and trim NAO data
load NAO_annual_DJFM.txt
NAO_annual_DJFM(NAO_annual_DJFM(:,1)<1955 | NAO_annual_DJFM(:,1)>2017,:)=[];

tmp=[1955:2017]';
weak_mwt=tmp>=1971 & tmp<=1975;
weaker_mwt=tmp>=2013 & tmp<=2017;
recent_mwt=time_final>1988;
tmp_time=1955.5:2017.5;
monthi=(1/6);
OCyr=time_final(temp_cont_final(1,:)<19);

f1=figure('units','centimeters','outerposition',[0.01 0.01 18 16]);
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
rectangle('position',[1978.5,min(pv_props_annual(:,1))+0.02,2,max(pv_props_annual(:,1))-...
    min(pv_props_annual(:,1))],'edgecolor','w','facecolor','w');
hold on
coeffs=polyfit(time_final,inpaint_nans(pv_data_final(:,1)),1);
fittedy=polyval(coeffs,1955:2017);
l1(1)=plot(1955:2017,fittedy,'--','color',rgb('black'),'linewidth',2);
coeffs=polyfit(time_final(recent_mwt),inpaint_nans(pv_data_final(recent_mwt,1)),1);
fittedy=polyval(coeffs,1988:2017);
l1(2)=plot(1988:2017,fittedy,':','color',rgb('black'),'linewidth',2);
plot(1955.5:2017.5,pv_props_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955.5:2017.5,pv_props_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--','ytick',16:21);
text(0.95,0.15,'d)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('\Theta (^oC)','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,pv_data_final(:,1)...
    ,2,[0.85 0.85 0.85],'filled');
ylim([15.8 max(pv_data_final(:,1))])
yl=ylim;
yl(1)=yl(1)+(max(yl)/1000);  yl(2)=yl(2)-(max(yl)/1000);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
uistack(s,'bottom');
% uistack(f3,'bottom')
uistack(f1,'bottom')
uistack(f2,'bottom')
uistack(l1,'top')

axes('position',[0.07 0.49 0.4 0.185]);
hold on
tmp_err=~isnan(pv_min_err(:,2));
jbfill(tmp_time(tmp_err),pv_min_err(tmp_err,2),flipud(pv_min_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),pv_min_err(tmp_err,2),flipud(pv_min_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978.5,min(pv_min_annual(:,1))+2e-12,2,max(pv_min_annual(:,1))-...
    min(pv_min_annual(:,1))],'edgecolor','w','facecolor','w');
hold on
coeffs=polyfit(time_final,inpaint_nans(low_pv_final(:,1)),1);
fittedy=polyval(coeffs,1955:2017);
l1(1)=plot(1955:2017,fittedy,'--','color',rgb('black'),'linewidth',2);
coeffs=polyfit(time_final(recent_mwt),inpaint_nans(low_pv_final(recent_mwt,1)),1);
fittedy=polyval(coeffs,1988:2017);
l1(2)=plot(1988:2017,fittedy,':','color',rgb('black'),'linewidth',2);
plot(1955.5:2017.5,pv_min_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955.5:2017.5,pv_min_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.15,'e)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('Core PV (m^{-1} s^{-1})','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final(low_pv_final(:,1)>0),low_pv_final(low_pv_final(:,1)>0,1)...
    ,2,[0.85 0.85 0.85],'filled');
ylim([-2e-12 max(low_pv_final(:,1))]);
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
uistack(l1,'top')

axes('position',[0.07 0.24 0.4 0.185]);
hold on
tmp_err=~isnan(intensity_err(:,2));
jbfill(tmp_time(tmp_err),intensity_err(tmp_err,2),flipud(intensity_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),intensity_err(tmp_err,2),flipud(intensity_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978.5,min(intensity_err(:,1))+2e-9,2,max(intensity_err(:,1))-...
    min(intensity_err(:,1))+2e-9],'edgecolor','w','facecolor','w');
hold on
coeffs=polyfit(time_final,inpaint_nans(worthington_intensity_final(:,1)),1);
fittedy=polyval(coeffs,1955:2017);
l1(1)=plot(1955:2017,fittedy,'--','color',rgb('black'),'linewidth',2);
coeffs=polyfit(time_final(recent_mwt),inpaint_nans(worthington_intensity_final(recent_mwt,1)),1);
fittedy=polyval(coeffs,1988:2017);
l1(2)=plot(1988:2017,fittedy,':','color',rgb('black'),'linewidth',2);
plot(1955.5:2017.5,intensity_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955.5:2017.5,intensity_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.9,'g)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('Intensity (s^{-1})','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,worthington_intensity_final(:,1)...
    ,2,[0.85 0.85 0.85],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
uistack(l1,'top')

axes('position',[0.07 0.74 0.4 0.185]);
hold on
tmp_err=~isnan(temp_thick_err(:,2));
jbfill(tmp_time(tmp_err),temp_thick_err(tmp_err,2),flipud(temp_thick_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),temp_thick_err(tmp_err,2),flipud(temp_thick_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978.5,min(temp_thick_err(:,1))+50,2,max(temp_thick_err(:,1))-...
    50],'edgecolor','w','facecolor','w');
hold on
coeffs=polyfit(time_final,inpaint_nans(stevens_thickness_final(:,1)),1);
fittedy=polyval(coeffs,1955:2017);
l1(1)=plot(1955:2017,fittedy,'--','color',rgb('black'),'linewidth',2);
coeffs=polyfit(time_final(recent_mwt),inpaint_nans(stevens_thickness_final(recent_mwt,1)),1);
fittedy=polyval(coeffs,1988:2017);
l1(2)=plot(1988:2017,fittedy,':','color',rgb('black'),'linewidth',2);
plot(1955.5:2017.5,temp_thick_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955.5:2017.5,temp_thick_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'gridlinestyle','--','ytick',0:100:500);
text(0.95,0.9,'c)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel('Thickness (m)','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,stevens_thickness_final(:,1)...
    ,2,[0.85 0.85 0.85],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
uistack(l1,'top')


axes('position',[0.55 0.49 0.4 0.185]);
hold on
tmp_err=~isnan(stmw_heat_content_err(:,2));
jbfill(tmp_time(tmp_err),stmw_heat_content_err(tmp_err,2),flipud(stmw_heat_content_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),stmw_heat_content_err(tmp_err,2),flipud(stmw_heat_content_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978.5,min(stmw_heat_content_err(:,1))+1e6,2,max(stmw_heat_content_err(:,1))-...
    min(stmw_heat_content_err(:,1))],'edgecolor','w','facecolor','w');
hold on
coeffs=polyfit(time_final,inpaint_nans(stmw_heat_content_final(:,1)),1);
fittedy=polyval(coeffs,1955:2017);
l1(1)=plot(1955:2017,fittedy,'--','color',rgb('black'),'linewidth',2);
coeffs=polyfit(time_final(recent_mwt),inpaint_nans(stmw_heat_content_final(recent_mwt,1)),1);
fittedy=polyval(coeffs,1988:2017);
l1(2)=plot(1988:2017,fittedy,':','color',rgb('black'),'linewidth',2);
plot(1955.5:2017.5,stmw_heat_content_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955.5:2017.5,stmw_heat_content_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.9,'f)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel({'STMW Integrated';'S.Enthalpy (J kg^{-1} m)'},'FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,stmw_heat_content_final(:,1)...
    ,2,[0.85 0.85 0.85],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/100);  yl(2)=yl(2)-(max(yl)/100);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
uistack(l1,'top')

axl=axes('position',[0.55 0.24 0.4 0.185]);
hold on
tmp_err=~isnan(heat_content_err(:,2));
jbfill(tmp_time(tmp_err),heat_content_err(tmp_err,2),flipud(heat_content_err(tmp_err,1)),'w','w',1,1);
jbfill(tmp_time(tmp_err),heat_content_err(tmp_err,2),flipud(heat_content_err(tmp_err,1)),...
    rgb('light red'),rgb('red'),1,0.5);
rectangle('position',[1978.5,min(heat_content_err(:,1))+5e4,2,max(heat_content_err(:,1))-...
    min(heat_content_err(:,1))],'edgecolor','w','facecolor','w');
hold on
coeffs=polyfit(time_final,inpaint_nans(heat_content_final(:,1)),1);
fittedy=polyval(coeffs,1955:2017);
l1(1)=plot(1955:2017,fittedy,'--','color',rgb('black'),'linewidth',2);
coeffs=polyfit(time_final(recent_mwt),inpaint_nans(heat_content_final(recent_mwt,1)),1);
fittedy=polyval(coeffs,1988:2017);
l1(2)=plot(1988:2017,fittedy,':','color',rgb('black'),'linewidth',2);
plot(1955.5:2017.5,heat_content_annual(:,1),'k','linewidth',1,'linejoin','round');
scatter(1955.5:2017.5,heat_content_annual(:,1),5,'w','filled','markeredgecolor','k');
axis tight
grid on
set(gca,'GridLineStyle','--');
text(0.95,0.15,'h)','Units','normalized','FontSize',8,'FontWeight','bold');
ylabel({'Surface Integrated';'S.Enthalpy (J kg^{-1} m)'},'FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8);
s=scatter(time_final,heat_content_final(:,1)...
    ,2,[0.85 0.85 0.85],'filled');
yl=ylim;
yl(1)=yl(1)+(max(yl)/500);  yl(2)=yl(2)-(max(yl)/500);
f1=patch([min(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) max(tmp_time(weak_mwt)) min(tmp_time(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) max(tmp_time(weaker_mwt)) min(tmp_time(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
box on
uistack(s,'bottom');
uistack(f1,'bottom')
uistack(f2,'bottom')
uistack(l1,'top')

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
xlim([1955 2018.5]);


%% STATS- Weak event
tmp=[1955:2017]';
weak_mwt=time_final>=1971 & time_final<=1976;
weaker_mwt=time_final>=2013 & time_final<=2018;

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

beep
pause

%% FIGURE- Contour
clear
lnes=lines;

% if ~exist('stevens_thickness_final','var')
load('20190504_MWdataset');
% end

pv_cont_final=pv_cont(:,(time>0));
% heat_content_final(heat_content_final>2.5e7)=NaN;
% heat_content_final=inpaint_nans(heat_content_final);
% aou_final(aou_final>100)=NaN;

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

% Create necessary datasets
% Make 2 monthly averaged temperature and limits
tmonth=zeros(120,113);
pmonth=zeros(120,113);
tgmonth=zeros(120,113);
wlimmonth=zeros(113,2);
elimmonth=zeros(113,2);
timecount=1989:1/6:max(time_final);

for i=1:length(timecount)
    ctime=timecount(i);
    tmonth(:,i)=nanmean(temp_cont_final(:,(time_final>=ctime & time_final<ctime+1/4)),2);
    pmonth(:,i)=nanmean(pv_cont_final(:,(time_final>=ctime & time_final<ctime+1/4)),2);
    tgmonth(:,i)=nanmean(tgrad_cont_final(:,(time_final>=ctime & time_final<ctime+1/4)),2);
%     wlimmonth(i,1)=nanmean(isotherm_limits_final((time_final>=ctime & time_final<ctime+1/4),1));
%     wlimmonth(i,2)=nanmean(isotherm_limits_final((time_final>=ctime & time_final<ctime+1/4),2));
end
% 
% timecount(isnan(wlimmonth(:,1)))=[];
% tmonth(:,(isnan(tmonth(100,:))))=[];
% tmonth=inpaint_nans(tmonth);
% pmonth(:,(isnan(pmonth(100,:))))=[];
% pmonth=inpaint_nans(pmonth);
% wlimmonth((isnan(wlimmonth(:,1))),:)=[];
% elimmonth((isnan(elimmonth(:,1))),:)=[];

% Annual datasets
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

for i=1988:1:2017
    % Stevens thickness
    tmp_data=stevens_thickness_final(time_final>=i & time_final<i+1);
    annuali=nanmean(tmp_data);
    temp_thick_annual=[temp_thick_annual;annuali i];
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    temp_thick_err=[temp_thick_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
    % Intensity
    tmp_data=worthington_intensity_final(time_final>i & time_final<i+1);
    intensi=nanmean(tmp_data);
    intensity_annual=[intensity_annual;intensi i];
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    intensity_err=[intensity_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
    % PV Min. 
    tmp_data=low_pv_final(time_final>i & time_final<i+1,1);
    var_annual=nanmean(tmp_data);
    pv_min_annual=[pv_min_annual;var_annual i];       
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    pv_min_err=[pv_min_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
    
    % PV min. properties
    tmp_data=pv_data_final(time_final>i & time_final<i+1,1);
    var_annual=nanmean(tmp_data);
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
    
    tmp_data=stmw_heat_content_final((time_final>i & time_final<i+1),:);
    var_annual=nanmean(tmp_data);
    stmw_heat_content_annual=[stmw_heat_content_annual;var_annual i]; 
    SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
    stmw_heat_content_err=[stmw_heat_content_err;nanmean(tmp_data)+ts*SEM]; % Confidence Intervals
end

% Load and define NAO data
load NAO_annual_DJFM.txt
NAO_annual_DJFM(NAO_annual_DJFM(:,1)<1988,:)=[];

% Isotherm
f1=figure('units','centimeters','outerposition',[0.01 0.01 12 10]);
ax1=axes('position',[.1 .6 .8 .35]);
contourf(timecount,5:5:600,tmonth,'LineColor','none','levelstep',0.1,'LevelList',[13.5:0.25:21.5]);
colormap(cmocean('Thermal'));
c=colorbar;
caxis([14 21])
ylabel(c,'\Theta (^oC)','FontSize',8,'FontWeight','bold');
set(gca,'ydir','reverse','FontSize',8);
hold on
[C,h]=contour(timecount,5:5:600,tmonth,[17 19],'color',...
    [0.8 0.8 0.8],'linewidth',2);
ylim([50 600]);
text(0.01,1.1,'a)','Units','normalized','FontSize',8,'FontWeight','bold')

% PV
ax1=axes('position',[.1 .1 .8 .35]);
contourf(timecount,5:5:600,pmonth,'LineColor','none','LevelList',...
    [1e-11:5e-12:2e-10]);
colormap(ax1,m_colmap('diverging'));
c=colorbar;
ylabel('Depth (m)','FontSize',8,'FontWeight','bold');
ylabel(c,'PV (m^{-1} s^{-1})','FontSize',8,'FontWeight','bold');
set(gca,'ydir','reverse','FontSize',8);
hold on
[C,h]=contour(timecount,5:5:600,pmonth,[1e-10 1e-10],'color',...
    [0.4 0.4 0.4],'linewidth',2);
% plot(time_final,movmean(low_pv_final(:,2),20),'k','linewidth',2);
ylim([50 600]);
text(0.01,1.1,'b)','Units','normalized','FontSize',8,'FontWeight','bold')

%% Figure- Formation and BVAL
clear
load 20190504_formation_dataset.mat
BVal=load('BV_trend_data.mat');

% convert OC volume to Svy if not done already
if OC_volume(1)>1000
    OC_volume=OC_volume/(1e6*60*60*24*365);
    vol_err=vol_err/(1e6*60*60*24*365);
end

f2=figure('units','centimeters','outerposition',[0.01 0.01 12 12]);
lnes=lines;
ax1=axes('position',[.11 .58 .8 .4]);
yyaxis left
jbfill(1989.5:2017.5,BATS_stmw.thick_err(:,2),flipud(BATS_stmw.thick_err(:,1)),...
    rgb('light gold'),rgb('gold'),1,0.5);
jbfill(1993.5:2018.5,thick_err(:,2),flipud(thick_err(:,1)),...
    rgb('light red'),rgb('red'),1,0.5);
hold on
plot(BATS_stmw.thick_annual(:,2)+0.5,BATS_stmw.thick_annual(:,1),'k-','linewidth',1.5,'linejoin','round');
scatter(BATS_stmw.thick_annual(:,2)+0.5,BATS_stmw.thick_annual(:,1),7,'w','filled','markeredgecolor','k');
% errorbar(BATS_stmw.intensity_annual(:,2)+0.5,BATS_stmw.intensity_annual(:,1),...
%     BATS_stmw.intensity_annual(:,1)-BATS_stmw.intensity_err(:,1),...
%     BATS_stmw.intensity_err(:,2)-BATS_stmw.intensity_annual(:,1),'color','k');
plot(years+0.5,sat_thick,'k--','linewidth',1.5,'linejoin','round');
scatter(years+0.5,sat_thick,7,'w','filled','markeredgecolor','k');
% errorbar(years+0.5,sat_intensity,sat_intensity(:,1)-OC_intensity_err(:,1),...
%     OC_intensity_err(:,2)-sat_intensity(:,1),'k');
set(gca,'Ycolor','k','fontsize',8);
ylabel({'Thickness (m)'},'Fontsize',8,'fontweight','bold');
grid on
yyaxis right
plot(years+0.5,OC_volume,'color',lnes(4,:),'linewidth',1.5,'linejoin','round');
scatter(years+0.5,OC_volume,7,'w','filled','markeredgecolor',lnes(4,:));
% errorbar(years+0.5,OC_volume,OC_volume-vol_err(:,1),vol_err(:,2)-OC_volume,...
%     'color',lnes(4,:));
set(gca,'Ycolor',lnes(4,:));
axis tight
grid on
set(gca,'GridLineStyle','--');
ylabel({'Volume (Svy)'},'Fontsize',8,'fontweight','bold');
hold on
coeffs=polyfit(years,OC_volume,1);
fittedx=polyval(coeffs,years);
plot(years,fittedx,'linestyle','-.','color',lnes(4,:),'linewidth',1.5);
text(0.9,0.9,'a)','Units','normalized','FontSize',8,'FontWeight','bold');
box on

ax2=axes('position',[.1 .08 .2 .4]);
plot(BVal.BVSthick_trend(:,5)/1000,BVal.BVSthick_trend(:,1),'k','linewidth',1.5);
axis tight
grid on
set(gca,'GridLineStyle','--');
hold on
set(gca,'fontsize',8,'xtick',0:250:1000);
ylabel('Thickness trend (m yr^{-1})','fontsize',8,'FontWeight','bold');
text(0.1,0.9,'b)','Units','normalized','FontSize',8,'FontWeight','bold');
plot([BVal.BVSthick_trend(2,5)/1000 BVal.BVSthick_trend(2,5)/1000],...
    [min(BVal.BVSthick_trend(:,1)) max(BVal.BVSthick_trend(:,1))],...
    'k--');

ax3=axes('position',[.4 .08 .2 .4]);
plot(BVal.BVintens_trend(:,5)/1000,BVal.BVintens_trend(:,1),'k','linewidth',1.5);
axis tight
grid on
set(gca,'GridLineStyle','--');
hold on
set(gca,'fontsize',8,'xtick',0:250:1000);
ylabel('Intensity trend (s^{-1} yr^{-1})','fontsize',8,'FontWeight','bold');
xlabel('Along-Transect Distance (km)','fontsize',8','FontWeight','bold');
text(0.1,0.9,'c)','Units','normalized','FontSize',8,'FontWeight','bold');
plot([BVal.BVintens_trend(2,5)/1000 BVal.BVintens_trend(2,5)/1000],...
    [min(BVal.BVintens_trend(:,1)) max(BVal.BVintens_trend(:,1))],...
    'k--');

ax4=axes('position',[.7 .08 .2 .4]);
plot(BVal.BVpv_trend(:,5)/1000,BVal.BVpv_trend(:,1),'k','linewidth',1.5);
axis tight
grid on
set(gca,'GridLineStyle','--');
hold on
set(gca,'fontsize',8,'xtick',0:250:1000);
ylabel({'Core PV trend';'(m^{-1} s^{-1} yr^{-1})'},'fontsize',8','FontWeight','bold');
text(0.8,0.9,'d)','Units','normalized','FontSize',8,'FontWeight','bold');
plot([BVal.BVpv_trend(2,5)/1000 BVal.BVpv_trend(2,5)/1000],...
    [min(BVal.BVpv_trend(:,1)) max(BVal.BVpv_trend(:,1))],...
    'k--');

for i=1:length(BVal.combo_log)
    
    axes(ax2);
    if BVal.combo_log(i,1)==1
        scatter(BVal.BVSthick_trend(i,5)/1000,BVal.BVSthick_trend(i,1),...
            20,rgb('neon blue'),'filled','markeredgecolor','k');
    else
        scatter(BVal.BVSthick_trend(i,5)/1000,BVal.BVSthick_trend(i,1),...
            20,'w','filled','markeredgecolor','k');
    end
    
    axes(ax3);
    if BVal.combo_log(i,2)==1
        scatter(BVal.BVintens_trend(i,5)/1000,BVal.BVintens_trend(i,1),...
            20,rgb('neon blue'),'filled','markeredgecolor','k');
    else
        scatter(BVal.BVintens_trend(i,5)/1000,BVal.BVintens_trend(i,1),...
            20,'w','filled','markeredgecolor','k');
    end
    
    axes(ax4);
    if BVal.combo_log(i,3)==1
        scatter(BVal.BVpv_trend(i,5)/1000,BVal.BVpv_trend(i,1),...
            20,rgb('neon blue'),'filled','markeredgecolor','k');
    else
       scatter(BVal.BVpv_trend(i,5)/1000,BVal.BVpv_trend(i,1),...
            20,'w','filled','markeredgecolor','k');
    end
%     drawnow
end

% Add BATS point onto plot
coeffs_thick=polyfit(BATS_stmw.time_final,...
    BATS_stmw.stevens_thickness_final,1);
coeffs_intens=polyfit(BATS_stmw.time_final,...
    BATS_stmw.worthington_intensity_final,1);
coeffs_pv=polyfit(BATS_stmw.time_final,...
    inpaint_nans(BATS_stmw.low_pv_final(:,1)),1);

scatter(ax2,BVal.BVpv_trend(2,5)/1000,coeffs_thick(1),...
            20,rgb('bright red'),'filled','markeredgecolor','k');
scatter(ax3,BVal.BVpv_trend(2,5)/1000,coeffs_intens(1),...
            20,rgb('bright red'),'filled','markeredgecolor','k');
scatter(ax4,BVal.BVpv_trend(2,5)/1000,coeffs_pv(1),...
            20,rgb('bright red'),'filled','markeredgecolor','k');        

%% STATS- Formation Stats
clc
coeffs=polyfit(years,OC_volume,1);
fittedx=polyval(coeffs,years);
[RHO, PVAL]=corr(BATS_stmw.intensity_annual(5:end,1),OC_volume(1:end-1));
disp(['BATS intensity-OC volume r=',num2str(RHO)]);
[RHO, PVAL]=corr(BATS_stmw.intensity_annual(5:end,1),sat_intensity(1:end-1));
disp(['BATS intensity-Model Intensity r=',num2str(RHO)]);

[RHO, PVAL]=corr(BATS_stmw.thick_annual(5:end,1),OC_volume(1:end-1));
disp(['BATS thickness-OC volume r=',num2str(RHO)]);
[RHO, PVAL]=corr(BATS_stmw.thick_annual(5:end,1),sat_thick(1:end-1));
disp(['BATS thicness-Model thickness r=',num2str(RHO)]);

[RHO, PVAL]=Mann_Kendall(OC_volume,0.05);
err=sqrt((sum((fittedx-OC_volume).^2))/(length(years-2)))/...
    sqrt(sum((years-nanmean(years)).^2));
change(1)=(((coeffs(1)-err)*length(years))/fittedx(1))*-100;
change(2)=(((coeffs(1)+err)*length(years))/fittedx(1))*-100;

% Check correlations
[r,pval]=corr(BATS_stmw.surf_heat_annual(5:end,1),OC_volume(1:end-1));
disp(['BATS surface heat-OC volume r=',num2str(r)]);

[r,pval]=corr(BATS_stmw.intensity_annual(5:end,1),OC_position(1:end-1));
disp(['BATS intensity-OC position r=',num2str(r)]);

[r,pval]=corr(OC_volume,OC_position);
disp(['OC volume-OC position r=',num2str(r)]);

% OC volume change
tmpa=sprintf('%2.2f',coeffs(1));
tmpb=sprintf('%2.1f',change(1));
tmpc=sprintf('%2.1f',change(2));
tmpd=sprintf('%2.2f',err);
disp([tmpa,' +- ',tmpd,' Svy year reduction in OC volume']);
disp([tmpc,'-',tmpb,'% reduction in OC volume']);

tmpa=sprintf('%2.1f',(length(years)*(-1*coeffs(1)+err)));
tmpb=sprintf('%2.1f',(length(years)*(-1*coeffs(1)-err)));
disp([tmpb,' - ',tmpa,' Svy overall reduction in OC volume']);

% tmpa=sprintf('%2.1f',(length(years)*(-1*coeffs(1)+err))/3.15e13);
% tmpb=sprintf('%2.1f',(length(years)*(-1*coeffs(1)-err))/3.15e13);
% disp([tmpb,' - ',tmpa,' Svy reduction in OC volume']);

coeffs=polyfit(years,sat_intensity,1);
fittedx=polyval(coeffs,years);
[RHO, PVAL]=corr(BATS_stmw.intensity_annual(5:end,1),sat_intensity(1:end-1));
[RHO, PVAL]=Mann_Kendall(sat_intensity,0.05);
% disp([num2str(((max(fittedx)-min(fittedx))/max(fittedx))*100),'% reduction in OC intensity']);

% coeffs=polyfit(years,formation_pv_min_temp(:,1),1);

% Check pentads of OC position
count=0;
OC_pos_pent=NaN(length(years-5),1);
for i=years(1):years(end)-5
    count=count+1;
    time_idx=years>i & years<i+5;
    OC_pos_pent(count)=nanmean(OC_position(time_idx));
end
[tmp, idx]=max(OC_pos_pent);
disp('The pentad with the northernmost 18^oC outcropping was:');
disp([num2str(years(idx)),'-',num2str(years(idx)+4),' at ',num2str(tmp),'^oN']);

[RHO, PVAL]=corr(BATS_stmw.intensity_annual(5:end,1),OC_position(1:end-1));
[RHO, PVAL]=corr(sat_intensity,OC_position);

[RHO, PVAL]=Mann_Kendall(OC_position,0.05);

beep
pause

%% TABLE- print statistics to file

HS=load('20190504_extended_dataset');
HS.pv_data_final(HS.pv_data_final(:,1)<2,1)=NaN;
fid=fopen('STMW_table.txt','w');
fprintf(fid,'%-31s%-31s%-41s%-41s%-41s%-31s%-31s%-31s%-31s\r\n',...
    'Parameter','Period','Mean and std. dev.','Change',...
    'Slope and std. error','n','r^2','p','%');

masks=[HS.time_final>1900 & HS.time_final<1988 HS.time_final>1988 & HS.time_final<2018 ...
    HS.time_final>2013 & HS.time_final<2018  HS.time_final>1900 & HS.time_final<2018 ...
    HS.time_final>1988 & HS.time_final<2013];
qs=struct('thick',HS.stevens_thickness_final,'intens',HS.worthington_intensity_final,...
    'pv',HS.low_pv_final(:,1));%,'temp',HS.pv_data_final(:,1),'surf_enth',HS.heat_content_final,'stmw_enth',HS.stmw_heat_content_final);
myfields=fieldnames(qs);

for i=1:size(masks,2)
    for ii=1:length(myfields)
        idx=logical(~isnan(qs.(myfields{ii})).*masks(:,i));
        tmp_prop=qs.(myfields{ii})(idx);
        tmp_time=HS.time_final(idx);
        
        coeffs=polyfit(tmp_time,tmp_prop,1);
        fittedx=linspace(min(tmp_time),max(tmp_time),sum(idx));
        fittedy=polyval(coeffs,fittedx);
        r=corr(fittedy',tmp_prop)^2;
        err=sqrt((sum((fittedy'-tmp_prop).^2))/(sum(idx)-2))/...
            sqrt(sum((tmp_time-nanmean(tmp_time)).^2));
        change(1)=(coeffs(1)-err)*length(unique(floor(tmp_time)));
        change(2)=(coeffs(1)+err)*length(unique(floor(tmp_time)));
        pct(1)=(((coeffs(1)+err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        pct(2)=(((coeffs(1)-err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        [~,sig]=Mann_Kendall(tmp_prop,0.05);
%         if i==3 && ii==1
%         keyboard
%         end
        fprintf(fid,'%-30s %-4.0f-%-25.0f %-10.4e+-%-30.4e %-15.4e-%-25.4e %-15.4e+-%-25.4e %-30.0f %-30.2f %-30.2f %-4.0f-%-25.0f\r\n',...
            myfields{ii},floor(min(tmp_time)),ceil(max(tmp_time)),mean(tmp_prop),...
            std(tmp_prop),change(1),change(2),coeffs(1),err,...
            sum(idx),r,sig,pct(1),pct(2));
    end
    
    fprintf(fid,'\r\n');
end
fclose(fid);
clc
fid=fopen('STMW_phys_table.txt','w');
fprintf(fid,'%-31s%-31s%-41s%-41s%-41s%-31s%-31s%-31s%-31s\r\n',...
    'Parameter','Period','Mean and std. dev.','Change',...
    'Slope and std. error','n','r^2','p','%');

masks=[HS.time_final>1900 & HS.time_final<1988 HS.time_final>1988 & HS.time_final<2018 ...
    HS.time_final>2013 & HS.time_final<2018  HS.time_final>1900 & HS.time_final<2018 ...
    HS.time_final>1988 & HS.time_final<2013];
qs=struct('temp',HS.pv_data_final(:,1),'STMWenth',HS.stmw_heat_content_final,...
    'Surfenth',HS.heat_content_final);%,'temp',HS.pv_data_final(:,1),'surf_enth',HS.heat_content_final,'stmw_enth',HS.stmw_heat_content_final);
myfields=fieldnames(qs);

for i=1:size(masks,2)
    for ii=1:length(myfields)
        idx=logical(~isnan(qs.(myfields{ii})).*masks(:,i));
        tmp_prop=qs.(myfields{ii})(idx);
        tmp_time=HS.time_final(idx);
        
        coeffs=polyfit(tmp_time,tmp_prop,1);
        fittedx=linspace(min(tmp_time),max(tmp_time),sum(idx));
        fittedy=polyval(coeffs,fittedx);
        r=corr(fittedy',tmp_prop)^2;
        err=sqrt((sum((fittedy'-tmp_prop).^2))/(sum(idx)-2))/...
            sqrt(sum((tmp_time-nanmean(tmp_time)).^2));
        change(1)=(coeffs(1)-err)*length(unique(floor(tmp_time)));
        change(2)=(coeffs(1)+err)*length(unique(floor(tmp_time)));
        pct(1)=(((coeffs(1)+err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        pct(2)=(((coeffs(1)-err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        [~,sig]=Mann_Kendall(tmp_prop,0.05);
%         if i==3 && ii==1
%         keyboard
%         end
        fprintf(fid,'%-30s %-4.0f-%-25.0f %-10.4e+-%-30.4e %-15.4e-%-25.4e %-15.4e+-%-25.4e %-30.0f %-30.2f %-30.2f %-4.0f-%-25.0f\r\n',...
            myfields{ii},floor(min(tmp_time)),ceil(max(tmp_time)),mean(tmp_prop),...
            std(tmp_prop),change(1),change(2),coeffs(1),err,...
            sum(idx),r,sig,pct(1),pct(2));
        
        if strcmp(myfields{ii},'temp') && sum(idx)>5000
            fprintf(['Per century temp. change:\n %3.2f - %3.2f \n'],...
                (coeffs(1)-err)*100,(coeffs(1)+err)*100);
             fprintf(['Mean STMW temp:\n %3.2f'],...
                nanmean(tmp_prop));
            fprintf(['Mean 2017 STMW temp:\n %3.2f'],...
                nanmean(tmp_prop(tmp_time>2017 & tmp_time<2018)));
        end
    end
    
    fprintf(fid,'\r\n');
end
fclose(fid);