%% BIG DATA STMW
addpath(genpath('C:\Users\samst\Documents'))
clear

data.bats=load('20191010_extended_dataset');
data.argo=load('argo_dataset');
data.EN4=load('EN4g10_BATS.mat');
data.ORA=load('ORAS4_BATS.mat');
structs=fieldnames(data);

load('EN4g10_all.mat','outcrop_grid');EN4_o=outcrop_grid;
load('ORAS4_all.mat','outcrop_grid');ORA_o=outcrop_grid;

% Clean BATS/HS data
data.bats.pv_cont_final=data.bats.pv_cont(:,(data.bats.time>0));
data.bats.low_pv_final(data.bats.low_pv_final(:,1)<0 |...
 data.bats.low_pv_final(:,1)>2e-10,1)=NaN;
data.bats.low_pv_final(:,1)=inpaint_nans(data.bats.low_pv_final(:,1));
data.bats.worthington_final=data.bats.worthington(data.bats.time>0);
data.bats.heat_content_final(data.bats.heat_content_final>2.5e7 |...
 data.bats.heat_content_final<9e6)=NaN;
data.bats.heat_content_final=inpaint_nans(data.bats.heat_content_final);
data.bats.stmw_heat_content_final(data.bats.heat_content_final<9e6)=NaN;
data.bats.stmw_heat_content_final=inpaint_nans(data.bats.stmw_heat_content_final);

data.bats.pv_data_final((data.bats.pv_data_final(:,1)==0),1)=NaN;
data.bats.pv_data_final(:,1)=inpaint_nans(data.bats.pv_data_final(:,1));
data.bats.pv_data_final((data.bats.pv_data_final(:,2)==0),2)=NaN;
data.bats.pv_data_final(:,2)=inpaint_nans(data.bats.pv_data_final(:,2));
data.bats.pv_data_final((data.bats.pv_data_final(:,3)==0),3)=NaN;
data.bats.pv_data_final(:,3)=inpaint_nans(data.bats.pv_data_final(:,3));
data.bats.pv_data_final((data.bats.pv_data_final(:,4)==0),4)=NaN;
data.bats.pv_data_final(:,4)=inpaint_nans(data.bats.pv_data_final(:,4));
data.bats.pv_data_final((data.bats.pv_data_final(:,5)==0),5)=NaN;
data.bats.pv_data_final(:,5)=inpaint_nans(data.bats.pv_data_final(:,5));

% Convert all time to same format
data.argo.time_final=cnv_date(data.argo.mtime);
data.EN4.time_final=data.EN4.month;
data.ORA.time_final=data.ORA.month;

data.argo.stevens_thickness_final=data.argo.stevens_thickness;
data.EN4.stevens_thickness_final=data.EN4.stevens_thickness;
data.ORA.stevens_thickness_final=data.ORA.stevens_thickness;

data.argo.worthington_intensity_final=data.argo.worthington_intensity;
data.EN4.worthington_intensity_final=data.EN4.worthington_intensity;
data.ORA.worthington_intensity_final=data.ORA.worthington_intensity;

data.argo.pv_data_final=data.argo.pv_data;
data.EN4.pv_data_final=data.EN4.pv_data;
data.ORA.pv_data_final=data.ORA.pv_data;

data.argo.low_pv_final=data.argo.low_pv;
data.EN4.low_pv_final=data.EN4.low_pv;
data.ORA.low_pv_final=data.ORA.low_pv;

%% Annual datasets
% if strcmp(question,'Y')
years=1955:2018;
% else
%     years=1955:2018;
% end
fldnom=fieldnames(data);
for ii=1:length(structs)
    
    data.(structs{ii}).temp_thick_annual=[];
    data.(structs{ii}).temp_thick_err=[];
    data.(structs{ii}).pv_min_annual=[];
    data.(structs{ii}).pv_min_err=[];
    data.(structs{ii}).pv_props_annual=[];
    data.(structs{ii}).pv_props_err=[];
    data.(structs{ii}).intensity_annual=[];
    data.(structs{ii}).intensity_err=[];
    
    for i=years(1):years(end)
                        
        % isotherm thickness
        tmp_data=data.(structs{ii}).stevens_thickness_final(data.(structs{ii}).time_final>=i &  data.(structs{ii}).time_final<i+1);
        annuali=nanmean(tmp_data);
        data.(structs{ii}).temp_thick_annual=[data.(structs{ii}).temp_thick_annual;annuali i];
        SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
        ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
        data.(structs{ii}).temp_thick_err=[data.(structs{ii}).temp_thick_err;annuali+ts*SEM]; % Confidence Intervals
        
        % Intensity
        tmp_data= data.(structs{ii}).worthington_intensity_final( data.(structs{ii}).time_final>i &  data.(structs{ii}).time_final<i+1);
        intensi=nanmean(tmp_data,1);
        data.(structs{ii}).intensity_annual=[data.(structs{ii}).intensity_annual;intensi i];
        SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
        ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
        data.(structs{ii}).intensity_err=[data.(structs{ii}).intensity_err;intensi+ts*SEM]; % Confidence Intervals
        
        % PV Min.
        tmp_data=data.(structs{ii}).low_pv_final(( data.(structs{ii}).time_final>i &  data.(structs{ii}).time_final<i+1),1);
        var_annual=nanmean(tmp_data,1);
        data.(structs{ii}).pv_min_annual=[data.(structs{ii}).pv_min_annual;var_annual i];
        SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
        ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
        data.(structs{ii}).pv_min_err=[data.(structs{ii}).pv_min_err;var_annual+ts*SEM]; % Confidence Intervals
        
        % PV min. properties
        tmp_data= data.(structs{ii}).pv_data_final(data.(structs{ii}).time_final>i &  data.(structs{ii}).time_final<i+1,1);
        var_annual=nanmean(tmp_data,1);
        data.(structs{ii}).pv_props_annual=[data.(structs{ii}).pv_props_annual;var_annual i];
        SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
        ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
        data.(structs{ii}).pv_props_err=[data.(structs{ii}).pv_props_err;var_annual+ts*SEM]; % Confidence Intervals
    end
end    

data.EN4.years=min(floor(data.EN4.time_final)):max(floor(data.EN4.time_final));
data.ORA.years=min(floor(data.ORA.time_final)):max(floor(data.ORA.time_final));

for i=1:length(data.EN4.years)
    EN4_outcrop_ts(i)=sum(sum(EN4_o(:,:,i),1),2)/(1e6*60*60*24*365);
end

for i=1:length(data.ORA.years)
    ORA_outcrop_ts(i)=sum(sum(ORA_o(:,:,i),1),2)/(1e6*60*60*24*365);
end

% Remove first few years(2001-2002) of argo data (lack of profiles)
idx=years<=2003;
data.argo.temp_thick_annual(idx,:)=NaN;data.argo.temp_thick_err(idx,:)=NaN;
data.argo.pv_min_annual(idx,:)=NaN;data.argo.pv_min_err(idx,:)=NaN;
data.argo.pv_props_annual(idx,:)=NaN;data.argo.pv_props_err(idx,:)=NaN;
data.argo.intensity_annual(idx,:)=NaN;data.argo.intensity_err(idx,:)=NaN;

    
%% Bin Argo profiles
%
% [N,Xedges,Yedges] = histcounts2(data.argo.lon,data.argo.lat,...
%     min(data.argo.lon):1:max(data.argo.lon),...
%     min(data.argo.lat):1:max(data.argo.lat));
% Xcent=Xedges(1:end-1)+0.5;
% Ycent=Yedges(1:end-1)+0.5;

x=min(data.argo.lon):1:max(data.argo.lon);
y=min(data.argo.lat):1:max(data.argo.lat);
[X,Y]=meshgrid(x,y);
N=zeros(size(X));

for i=1:length(x)
    for ii=1:length(y)
        N(ii,i)=sum(data.argo.lon>X(ii,i)-0.5 & data.argo.lon<X(ii,i)+0.5 &...
            data.argo.lat>Y(ii,i)-0.5 & data.argo.lat<Y(ii,i)+0.5);
    end
end



%% PLOT
tmp=1955:2018;

weak_mwt=tmp>=1971 & tmp<=1975;
weaker_mwt=tmp>=2014 & tmp<=2018;

Fgm=figure('units','centimeters','outerposition',[0 0 9 16]);
set(gcf,'color','w');

ax(1)=subplot(5,1,1);
hold on
ax(1).XAxis.Visible='off';

ax(2)=subplot(5,1,2);
hold on
ax(2).XAxis.Visible='off';

ax(3)=subplot(5,1,3);
hold on
ax(3).XAxis.Visible='off';

ax(4)=subplot(5,1,4);
hold on
ax(4).XAxis.Visible='off';

ax(5)=subplot(5,1,5);
hold on

cm=linspecer(length(structs));
years=1955:2018;

for ii=1:length(structs)
    
    axes(ax(1));
    idx=~isnan(data.(structs{ii}).temp_thick_err(:,2));
%     jbfill(years(idx),data.(structs{ii}).temp_thick_err(idx,2),...
%         flipud(data.(structs{ii}).temp_thick_err(idx,1)),'w','w',1,1);
    jbfill(years(idx),data.(structs{ii}).temp_thick_err(idx,2),...
        flipud(data.(structs{ii}).temp_thick_err(idx,1)),...
        cm(ii,:),cm(ii,:),1,0.2);
    hold on
    
%     errorbar(years(idx),data.(structs{ii}).temp_thick_annual(:,1),...
%         data.(structs{ii}).temp_thick_err(:,1).*-1,...
%         data.(structs{ii}).temp_thick_err(:,2),...
%         'o','MarkerSize',2,'MarkerEdgeColor',cm(ii,:),'MarkerFaceColor','w');
     
    axes(ax(2));
    idx=~isnan(data.(structs{ii}).pv_min_err(:,2));
%     jbfill(years(idx),data.(structs{ii}).pv_min_err(idx,2),...
%         flipud(data.(structs{ii}).pv_min_err(idx,1)),'w','w',1,1);
    jbfill(years(idx),data.(structs{ii}).pv_min_err(idx,2),...
        flipud(data.(structs{ii}).pv_min_err(idx,1)),...
        cm(ii,:),cm(ii,:),1,0.2);
    hold on
    
    axes(ax(3));
    idx=~isnan(data.(structs{ii}).intensity_err(:,2));
%     jbfill(years(idx),data.(structs{ii}).intensity_err(idx,2),...
%         flipud(data.(structs{ii}).intensity_err(idx,1)),'w','w',1,1);
    jbfill(years(idx),data.(structs{ii}).intensity_err(idx,2),...
        flipud(data.(structs{ii}).intensity_err(idx,1)),...
        cm(ii,:),cm(ii,:),1,0.2);
    hold on
    
    axes(ax(4));
    idx=~isnan(data.(structs{ii}).pv_props_err(:,2));
%     jbfill(years(idx),data.(structs{ii}).pv_props_err(idx,2),...
%         flipud(data.(structs{ii}).pv_props_err(idx,1)),'w','w',1,1);
    jbfill(years(idx),data.(structs{ii}).pv_props_err(idx,2),...
        flipud(data.(structs{ii}).pv_props_err(idx,1)),...
        cm(ii,:),cm(ii,:),1,0.2);
    hold on
end


for ii=1:length(structs)
    
    axes(ax(1));
    l(ii)=plot(1955:2018,movmean(...
        data.(structs{ii}).temp_thick_annual(:,1),3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(ii,:));
         
%     errorbar(1955:2018,data.(structs{ii}).temp_thick_annual(:,1),...
%         data.(structs{ii}).temp_thick_err(:,1).*-1,...
%         data.(structs{ii}).temp_thick_err(:,2),...
%         'o','MarkerSize',2,'MarkerEdgeColor',cm(ii,:),'MarkerFaceColor','w');
     
    axes(ax(2));
    plot(1955:2018,movmean(...
        data.(structs{ii}).pv_min_annual(:,1),3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(ii,:));
    
    axes(ax(3));
    plot(1955:2018,movmean(...
        data.(structs{ii}).intensity_annual(:,1),3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(ii,:));
    
    axes(ax(4));
    plot(1955:2018,movmean(...
        data.(structs{ii}).pv_props_annual(:,1),3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(ii,:));
end

axes(ax(5));
plot(data.EN4.years,movmean(...
        EN4_outcrop_ts,3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(3,:));
 plot(data.ORA.years,movmean(...
        ORA_outcrop_ts,3,'omitnan'),...
        'linewidth',1.5,'linejoin','round','color',cm(4,:));   

axes(ax(1));
axis tight
yl=ylim;
ylim([0 max(yl)]);
ax(1).Color='none';
ylabel('Thickness (m)','FontSize',8,'FontWeight','bold');
ax(1).YGrid='on';
ax(1).GridAlpha = 0.06;
xtick([0 100 200]);
text(.04,.12,'a)','Units','normalized','FontSize',8,'FontWeight','bold')

axes(ax(2));
axis tight
ax(2).Color='none';
ylabel({'Core PV (s^{-1} m^{-1})'},'FontSize',8,'FontWeight','bold');
ax(2).YGrid='on';
ax(2).GridAlpha = 0.06;
text(.04,.12,'b)','Units','normalized','FontSize',8,'FontWeight','bold')

axes(ax(3));
axis tight
ax(3).Color='none';
ylabel({'Intensity';'(s^{-1})'},'FontSize',8,'FontWeight','bold');
ax(3).YGrid='on';
ax(3).GridAlpha = 0.06;
text(.04,.12,'c)','Units','normalized','FontSize',8,'FontWeight','bold')

axes(ax(4));
axis tight
ax(4).Color='none';
ylabel('\Theta (^oC)','FontSize',8,'FontWeight','bold');
ax(4).YGrid='on';
ax(4).GridAlpha = 0.06;
text(.04,.12,'d)','Units','normalized','FontSize',8,'FontWeight','bold')

axes(ax(5));
axis tight
ax(5).Color='none';
ylabel({'Outcropping';'Volume (Svy)'},'FontSize',8,'FontWeight','bold');
ax(5).YGrid='on';
ax(5).GridAlpha = 0.06;
text(.04,.12,'e)','Units','normalized','FontSize',8,'FontWeight','bold')

iax=axes;
plot(1955:2018,NaN(length(1955:2018),1));
yl=ylim;
yl(1)=yl(1)+(max(yl)/1000);  yl(2)=yl(2)-(max(yl)/1000);
f1=patch([min(tmp(weak_mwt)) max(tmp(weak_mwt)) max(tmp(weak_mwt)) min(tmp(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp(weaker_mwt)) max(tmp(weaker_mwt)) max(tmp(weaker_mwt)) min(tmp(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
linkaxes([iax ax],'x');
iax.XAxis.Visible='off';
iax.YAxis.Visible='off';
iax.XGrid='on';
iax.XMinorGrid='on';
iax.MinorGridLineStyle='-';
iax.GridAlpha = 0.075;
iax.MinorGridAlpha = 0.04;
% y=legend(l,{'BATS/HS','Argo','EN4g10','ORAS4'},'FontSize',6,'Location','best');
% title(y,'Dataset');
uistack(iax,'bottom');

%% When legend is adjusted...

export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\all_ts.eps -deps -painters
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\all_ts.pdf -dpdf -painters

% close all

%% Stats for table 2

tmp=[1955:2018]';

for i=1:length(structs)
    
    weak_mwt=data.(structs{i}).time_final>=1971 & data.(structs{i}).time_final<=1976;
    weaker_mwt=data.(structs{i}).time_final>=2014 & data.(structs{i}).time_final<2019;
    
    disp(structs(i));
    tmp_data=data.(structs{i}).stevens_thickness_final(:,1);
    SEM=nanstd(tmp_data(weak_mwt))/sqrt(length(tmp_data(weak_mwt))); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data(weak_mwt))-1); % T-Score
    quick_err(1,:)=nanmean(tmp_data(weak_mwt))+ts*SEM; % Confidence Intervals
    disp(['Weak 70s event had ',num2str(quick_err(1,1)),'-',num2str(quick_err(1,2)),...
        ' m mean thickness, n=' num2str(sum(weak_mwt))]);
    
    SEM=nanstd(tmp_data(weaker_mwt))/sqrt(length(tmp_data(weaker_mwt))); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data(weaker_mwt))-1); % T-Score
    quick_err(2,:)=nanmean(tmp_data(weaker_mwt))+ts*SEM; % Confidence Intervals
    disp(['Weak recent event had ',num2str(quick_err(2,1)),'-',num2str(quick_err(2,2)),...
        ' m mean thickness, n=' num2str(sum(weaker_mwt))]);
    
    tmp_data=data.(structs{i}).worthington_intensity_final(:,1);
    SEM=nanstd(tmp_data(weak_mwt))/sqrt(length(tmp_data(weak_mwt))); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data(weak_mwt))-1); % T-Score
    quick_err(3,:)=nanmean(tmp_data(weak_mwt))+ts*SEM; % Confidence Intervals
    disp(['Weak 70s event had ',num2str(quick_err(3,1)),'-',num2str(quick_err(3,2)),...
        ' s-1 mean intensity, n=' num2str(sum(weak_mwt))]);
    
     SEM=nanstd(tmp_data(weaker_mwt))/sqrt(length(tmp_data(weaker_mwt))); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data(weaker_mwt))-1); % T-Score
    quick_err(4,:)=nanmean(tmp_data(weaker_mwt))+ts*SEM; % Confidence Intervals
    disp(['Weak recent event had ',num2str(quick_err(4,1)),'-',num2str(quick_err(4,2)),...
        ' s-1 mean intensity, n=' num2str(sum(weaker_mwt))]);
    
    tmp_data=data.(structs{i}).low_pv_final(:,1);
    SEM=nanstd(tmp_data(weak_mwt))/sqrt(length(tmp_data(weak_mwt))); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data(weak_mwt))-1); % T-Score
    quick_err(5,:)=nanmean(tmp_data(weak_mwt))+ts*SEM; % Confidence Intervals
    disp(['Weak 70s event had ',num2str(quick_err(5,1)),'-',num2str(quick_err(5,2)),...
        ' m-1 s-1 mean PV, n=' num2str(sum(weak_mwt))]);
    
     SEM=nanstd(tmp_data(weaker_mwt))/sqrt(length(tmp_data(weaker_mwt))); % Standard Error
    ts=tinv([0.025  0.975],length(tmp_data(weaker_mwt))-1); % T-Score
    quick_err(6,:)=nanmean(tmp_data(weaker_mwt))+ts*SEM; % Confidence Intervals
    disp(['Weak recent event had ',num2str(quick_err(6,1)),'-',num2str(quick_err(6,2)),...
        ' m-1 s-1 mean PV, n=' num2str(sum(weaker_mwt))]);
    
end



%% Map
fGMs=figure;
set(gcf,'color','w');
lon_lim=[-76 -44];
lat_lim=[29 43];
m_proj('lambert','lon',lon_lim,'lat',lat_lim,'rect','on');   % Projection
hold on
RGB=rgb('light grey');
% Contour formation zones
[C,h]=m_contourf(X,Y,N,'linestyle','none');
m_scatter(data.argo.lon,data.argo.lat,1.5,'k','filled');
cm=cmocean('amp');cm(1,:)=[1,1,1];
colormap(fGMs,cm);
[ax,h]=m_contfbar(-0.04,[0.2 0.8],C,h,'axfrac',.02,'fontsize',8,'fontweight','bold');
ylabel(ax,{'Profile';'Count'},'rotation',0,'color','k','fontsize',7);
m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',7,'box','on');
m_scatter(-64.16666,31.666,20,'k','filled');
m_text(-64.16666+0.5,31.666-0.1,'BATS','color','k','fontsize',6);
% m_scatter(-64.5,32.16666,25,'w','filled');
% m_text(-64.5+0.5,32.16666+0.1,'Hydrostation ''S''','color','k','fontweight','bold','fontsize',8);
%%
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Argo_map.pdf -dpdf -painters
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Argo_map.tif -dtif
