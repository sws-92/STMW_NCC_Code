%% Map reanalysis data
clear
addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))

ORA=load('ORAS4_all.mat');
tmp=load('ORAS4_data.mat','time');ORA.all_time=tmp.time;
EN4=load('EN4g10_all.mat');
load 20191028_BVstations

ORA.time=1958:.25:2018;
EN4.time=1955:.25:2019;
ORA.years=unique(floor(ORA.time));
EN4.years=unique(floor(EN4.time));

% Create annual values
ORA.annual_thick=NaN([size(ORA.lon_grid) length(ORA.years)]);
ORA.annual_intens=NaN([size(ORA.lon_grid) length(ORA.years)]);
ORA.annual_pv=NaN([size(ORA.lon_grid) length(ORA.years)]);
ORA.annual_stemp=NaN([size(ORA.lon_grid) length(ORA.years)]);

EN4.annual_thick=NaN([size(EN4.lon_grid) length(EN4.years)]);
EN4.annual_intens=NaN([size(EN4.lon_grid) length(EN4.years)]);
EN4.annual_pv=NaN([size(EN4.lon_grid) length(EN4.years)]);
EN4.annual_stemp=NaN([size(EN4.lon_grid) length(EN4.years)]);

for i=1:length(EN4.years)
    idx=EN4.time>EN4.years(i) & EN4.time<EN4.years(i)+1;
    EN4.annual_thick(:,:,i)=mean(EN4.thick_grid(:,:,idx),3,'omitnan');
    EN4.annual_intens(:,:,i)=mean(EN4.intensity_grid(:,:,idx),3,'omitnan');
    EN4.annual_pv(:,:,i)=mean(EN4.PV_grid(:,:,idx),3,'omitnan');
    EN4.annual_stemp(:,:,i)=mean(EN4.stemp_grid(:,:,idx),3,'omitnan');
end

for i=1:length(ORA.years)
    idx=ORA.time>ORA.years(i) & ORA.time<ORA.years(i)+1;
    ORA.annual_thick(:,:,i)=mean(ORA.thick_grid(:,:,idx),3,'omitnan');
    ORA.annual_intens(:,:,i)=mean(ORA.intensity_grid(:,:,idx),3,'omitnan');
    ORA.annual_pv(:,:,i)=mean(ORA.PV_grid(:,:,idx),3,'omitnan');
    ORA.annual_stemp(:,:,i)=mean(ORA.stemp_grid(:,:,idx),3,'omitnan');
end

load 20191028_BVstations
[x,y,z]=meshgrid(ORA.lon_grid(1,:),ORA.lat_grid(:,1),ORA.years);
[x2,y2,z2]=meshgrid(EN4.lon_grid(1,:),EN4.lat_grid(:,1),EN4.years);
[xn,yn,zn]=meshgrid(stations.coords_list(:,1),stations.coords_list(:,2),ORA.years);
[xn2,yn2,zn2]=meshgrid(stations.coords_list(:,1),stations.coords_list(:,2),EN4.years);

ORA_int=interp3(x,y,z,ORA.annual_thick,xn,yn,zn);
EN4_int=interp3(x2,y2,z2,EN4.annual_thick,xn2,yn2,zn2);

for i=1:length(stations.coords_list(:,1))
    stations.(['ORA_' num2str(i)])=squeeze(ORA_int(i,i,:));
    stations.(['EN4_' num2str(i)])=squeeze(EN4_int(i,i,:));
end

for i=1:length(EN4.years)
    EN4.outcrop_ts(i)=sum(sum(EN4.outcrop_grid(:,:,i),1),2);
end

for i=1:length(ORA.years)
    ORA.outcrop_ts(i)=sum(sum(ORA.outcrop_grid(:,:,i),1),2);
end
    
BATS=load('20191010_extended_dataset.mat','time_final','stevens_thickness_final');

BATS.temp_thick_annual=[];
BATS.temp_thick_err=[];
for i=floor(BATS.time_final(1)):floor(BATS.time_final(end))-1
        % isotherm thickness
        tmp_data=BATS.stevens_thickness_final(BATS.time_final>=i &  BATS.time_final<i+1);
        annuali=nanmean(tmp_data);
        BATS.temp_thick_annual=[BATS.temp_thick_annual;annuali i];
        SEM=nanstd(tmp_data)/sqrt(length(tmp_data)); % Standard Error
        ts=tinv([0.025  0.975],length(tmp_data)-1); % T-Score
        BATS.temp_thick_err=[BATS.temp_thick_err;ts*SEM]; % Confidence Intervals
end
        
    
%% Figure 1-ORA
all_temp=mean(ORA.wtemp_grid,3,'omitnan');
rec_temp=mean(ORA.wtemp_grid(:,:,ORA.years>=2013 & ORA.years<=2017),3,'omitnan');

% Set up figure
f1=figure('units','centimeters','outerposition',[0.01 0.01 18 13]);
% ax1=axes('position',[0.5 0 1 0.5]);
colormap(cmocean('Thermal'));
lon_lim=[-80 -40];
lat_lim=[19.5 45];
m_proj('lambert','lon',lon_lim,'lat',lat_lim,'rect','on');   % Projection
hold on
RGB=rgb('light grey');

% % Find formation box idxs
% lon_idx=lon_grid>-75 & lon_grid<-45;
% lat_idx=lat_grid>30 & lat_grid<41.5;

% Contour formation zones
[C,h]=m_contourf(ORA.lon_grid,ORA.lat_grid,all_temp,'LineColor','w',...
    'levelstep',1,'LevelList',2:1:23,'LineColor','none');
[ax,h]=m_contfbar(0.05,[0.2 0.8],C,h,'axfrac',.02,'FontSize',8,'fontweight','bold');
ylabel(ax,{'SST';'(^oC)'},'rotation',0,'color','k');
[C,c]=m_contour(ORA.lon_grid,ORA.lat_grid,all_temp,[17 19],'w--','linewidth',1.5);
[B,b]=m_contour(ORA.lon_grid,ORA.lat_grid,all_temp,[18 18],'w-.','linewidth',1.5);
[F,f]=m_contour(ORA.lon_grid,ORA.lat_grid,rec_temp,[17 19],'k--','linewidth',1);
[F,f]=m_contour(ORA.lon_grid,ORA.lat_grid,rec_temp,[18 18],'k-.','linewidth',1);

m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','FontSize',8,'box','on');
m_line([ORA.lon_grid(1,:) ORA.lon_grid(:,end)' fliplr(ORA.lon_grid(1,:)) ORA.lon_grid(:,1)'],...
    [ORA.lat_grid(end,:) fliplr(ORA.lat_grid(:,end)') ORA.lat_grid(1,:) ORA.lat_grid(:,end)'],...
    'color','k','linestyle','-.','linewidth',1.5);
[x,y]=meshgrid(-75:-45,30:41);
m_line([x(1,:) x(:,end)' fliplr(x(1,:)) x(:,1)'],...
    [y(end,:) fliplr(y(:,end)') y(1,:) y(:,end)'],...
    'color','k','linestyle','--','linewidth',1.5);
m_scatter(-64.16666,31.666,25,'k','filled');
m_text(-64.16666+0.5,31.666-0.1,'BATS','color','k','fontweight','bold','Fontsize',7);
m_scatter(-64.5,32.16666,25,'k','filled');
m_text(-64.5+0.5,32.16666+0.1,'Hydrostation ''S''','color','k','fontweight','bold','Fontsize',7);
% text(0.9,0.9,'a)','Units','normalized','FontSize',8,'FontWeight','bold');

m_scatter(stations.coords_list(:,1),stations.coords_list(:,2),15,...
    'k','filled');
m_text(stations.coords_list([1 3:7],1)+0.3,stations.coords_list([1 3:7],2),...
    {'BV1','BV3','BV4','BV5','BV6','BV7'},...
    'color','k','Fontsize',7);
%%
set(gcf,'color','w');
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure1_EN4S4.pdf -dpdf -painters

%% Figure 1- EN4
all_temp=mean(EN4.wtemp_grid(:,:,EN4.years<2013),3,'omitnan');
rec_temp=mean(EN4.wtemp_grid(:,:,EN4.years>=2014 & EN4.years<=2018),3,'omitnan');

% Set up figure
f1=figure('units','centimeters','outerposition',[0.01 0.01 18 13],'color','w');
% ax1=axes('position',[0.5 0 1 0.5]);
colormap(cmocean('Thermal'));
lon_lim=[-80 -40];
lat_lim=[19.5 45];
m_proj('lambert','lon',lon_lim,'lat',lat_lim,'rect','on');   % Projection
hold on
RGB=rgb('light grey');

% % Find formation box idxs
% lon_idx=lon_grid>-75 & lon_grid<-45;
% lat_idx=lat_grid>30 & lat_grid<41.5;

% Contour formation zones
[C,h]=m_contourf(EN4.lon_grid,EN4.lat_grid,all_temp,'LineColor','w',...
    'levelstep',1,'LevelList',2:1:23,'LineColor','none');
[ax,h]=m_contfbar(-0.015,[0.1 0.9],C,h,'axfrac',.03,'FontSize',8,'fontweight','bold');
ylabel(ax,{'SST';'(^oC)'},'rotation',0,'color','k');
[C,c]=m_contour(EN4.lon_grid,EN4.lat_grid,all_temp,[17 19],'w--','linewidth',1.5);
[B,b]=m_contour(EN4.lon_grid,EN4.lat_grid,all_temp,[18 18],'w-.','linewidth',1.5);
[F,f]=m_contour(EN4.lon_grid,EN4.lat_grid,rec_temp,[17 19],'k--','linewidth',1);
[F,f]=m_contour(EN4.lon_grid,EN4.lat_grid,rec_temp,[18 18],'k-.','linewidth',1);

m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','FontSize',8,'box','on');
m_line([EN4.lon_grid(1,:) EN4.lon_grid(:,end)' fliplr(EN4.lon_grid(1,:)) EN4.lon_grid(:,1)'],...
    [EN4.lat_grid(end,:) fliplr(EN4.lat_grid(:,end)') EN4.lat_grid(1,:) EN4.lat_grid(:,end)'],...
    'color','k','linestyle','-.','linewidth',1.5);
[x,y]=meshgrid(-75:-45,30:41);
m_line([x(1,:) x(:,end)' fliplr(x(1,:)) x(:,1)'],...
    [y(end,:) fliplr(y(:,end)') y(1,:) y(:,end)'],...
    'color','k','linestyle','--','linewidth',1.5);
m_scatter(-64.16666,31.666,25,'k','filled');
m_text(-64.16666+0.5,31.666-0.1,'BATS','color','k','fontweight','bold','Fontsize',7);
m_scatter(-64.5,32.16666,25,'k','filled');
m_text(-64.5+0.5,32.16666+0.1,'Hydrostation ''S''','color','k','fontweight','bold','Fontsize',7);
% text(0.9,0.9,'a)','Units','normalized','FontSize',8,'FontWeight','bold');

m_scatter(stations.coords_list(:,1),stations.coords_list(:,2),15,...
    'k','filled');
m_text(stations.coords_list([1 3:8],1)+0.3,stations.coords_list([1 3:8],2),...
    {'BV1','BV3','BV4','BV5','BV6','BV7','BV8'},...
    'color','k','Fontsize',7);

%%
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure1_EN4.pdf -dpdf -painters
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\EN4_map.eps -deps

%% Figure 6 - ORA 

f1=figure('units','centimeters','outerposition',[0.01 0.01 18 18]);
set(gcf,'color','w');
ax1=axes('position',[-0.05 0.45 0.8 0.45]);
lon_lim=[-80 -40];
lat_lim=[19.5 45];
m_proj('lambert','lon',lon_lim,'lat',lat_lim,'rect','on');   % Projection
hold on
RGB=rgb('light grey');

ORA.thick_change=mean(ORA.annual_thick(:,:,ORA.years<2014),3,'omitnan')-...
    mean(ORA.annual_thick(:,:,ORA.years>=2014),3,'omitnan');

% Contour formation zones
[C,h]=m_contourf(ORA.lon_grid,ORA.lat_grid,ORA.thick_change,0:25:250,'LineColor','none');
cm=cmocean('amp');cm(1,:)=[1,1,1];%;1,1,1;1,1,1];
colormap(ax1,cm);
[ax,h]=m_contfbar(ax1,[0.3 0.86],1.05,C,h,'axfrac',.03,'FontSize',8,...
    'fontweight','bold','endpiece','no');
text(0.5,1,{'Thickness';'change (m)'},'rotation',0,'color','k','FontSize',8,...
    'units','normalized','fontweight','bold');

m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','FontSize',8,'box','on',...
    'linest',':','color',[0.7,0.7,0.7]);
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','FontSize',8,'box','on');

m_line([ORA.lon_grid(1,:) ORA.lon_grid(:,end)' fliplr(ORA.lon_grid(1,:)) ORA.lon_grid(:,1)'],...
    [ORA.lat_grid(end,:) fliplr(ORA.lat_grid(:,end)') ORA.lat_grid(1,:) ORA.lat_grid(:,end)'],...
    'color','k','linestyle','-.','linewidth',1.5);
[x,y]=meshgrid(-75:-45,30:41);
m_line([x(1,:) x(:,end)' fliplr(x(1,:)) x(:,1)'],...
    [y(end,:) fliplr(y(:,end)') y(1,:) y(:,end)'],...
    'color','k','linestyle','--','linewidth',1.5);

cm=linspecer(length(stations.coords_list(:,1))+1);
cm(1,:)=[];
for i=1:length(stations.coords_list(:,1))
    m_scatter(stations.coords_list(i,1),stations.coords_list(i,2),30,...
        'filled','markerfacecolor',cm(i,:),'markeredgecolor','k');
end

axpos=linspace(0.875,0.03,length(stations.coords_list(:,1)));

for i=1:length(axpos)-1
    ax2(i)=axes('position',[0.675 axpos(i) 0.3 0.1]);
    hold on
    ax2(i).XAxis.Visible='off';
    ax2(i).YAxis.FontSize=6;
    ax2(i).Color='none';
end
ax2(i+1)=axes('position',[0.675 axpos(end) 0.3 0.1]);
hold on
ax2(i+1).XAxis.FontSize=6;
ax2(i+1).YAxis.FontSize=6;
ax2(i+1).Color='none';


for i=1:length(stations.coords_list(:,1))
    axes(ax2(i));
    tmp_thick=stations.(['Sthick' num2str(i)]);
    tmp_ORA=stations.(['ORA_' num2str(i)]);
    tmp_EN4=stations.(['EN4_' num2str(i)]);
    
    
    O_l(i)=plot(ax2(i),ORA.years,tmp_ORA,'color',rgb('light red'),...
        'linewidth',2);
    E_l(i)=plot(ax2(i),EN4.years,tmp_EN4,'color',rgb('light blue'),...
        'linewidth',2);
    s(i)=scatter(ax2(i),tmp_thick(:,2),tmp_thick(:,1),15,'filled',...
        'markerfacecolor',cm(i,:),'markeredgecolor','k','markerfacealpha',0.5);
    axis tight
    ylim([0 400]);
    ax2(i).YGrid='on';
    ax2(i).GridAlpha = 0.1;
end    

tmp=1955:2018;
weak_mwt=tmp>=1971 & tmp<=1975;
weaker_mwt=tmp>=2014 & tmp<=2018;

iax=axes('position',[0.675 axpos(end) 0.3 axpos(1)+0.075]);
plot(1955:2018,NaN(length(1955:2018),1));
yl=ylim;
yl(1)=yl(1)+(max(yl)/1000);  yl(2)=yl(2)-(max(yl)/1000);
f1=patch([min(tmp(weak_mwt)) max(tmp(weak_mwt)) max(tmp(weak_mwt)) min(tmp(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp(weaker_mwt)) max(tmp(weaker_mwt)) max(tmp(weaker_mwt)) min(tmp(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
iax.XAxis.Visible='off';
iax.YAxis.Visible='off';
iax.XGrid='on';
iax.XMinorGrid='on';
iax.MinorGridLineStyle='-';
iax.GridAlpha = 0.075;
iax.MinorGridAlpha = 0.04;
% y=legend(l,{'BATS/HS','Argo','EN4g10','ORAS4'},'Fontsize',7,'Location','best');
% title(y,'Dataset');
uistack(iax,'bottom');

ax3=axes('position',[0.072 0.1 0.525 0.275]);
box on
hold on
set(gca,'fontsize',7);
plot(ORA.years,smooth(ORA.outcrop_ts,3),'linewidth',2,...
    'color',rgb('light red'));
plot(EN4.years,smooth(EN4.outcrop_ts,3),'linewidth',2,...
    'color',rgb('light blue'));
yyaxis right
plot(BATS.temp_thick_annual(:,2),smooth(BATS.temp_thick_annual(:,1),3),':',...
    'linewidth',1,'color',rgb('black'));
axis tight
set(gca,'Ycolor','k');
ylabel('Thickness (m)','fontsize',8);

yyaxis left
axis tight
yl=ylim;ylim([0 max(ylim)]);
% xlim([min
grid on
ax3.GridAlpha = 0.1;
linkaxes([iax ax2 ax3],'x');
ylabel('Outcropping Volume (Svy)','fontsize',8);

%%
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure5_ORA.pdf -dpdf -painters
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure5_ORA.eps -deps -painters

%% Figure 6 - EN4
f=figure('units','centimeters','outerposition',[0.01 0.01 18 16]);
set(gcf,'color','w');
ax1=axes('position',[-0.05 0.45 0.8 0.45]);
lon_lim=[-80 -40];
lat_lim=[19.5 45];
m_proj('lambert','lon',lon_lim,'lat',lat_lim,'rect','on');   % Projection
hold on
RGB=rgb('light grey');

EN4.thick_change=mean(EN4.annual_thick(:,:,EN4.years<2014),3,'omitnan')-...
    mean(EN4.annual_thick(:,:,EN4.years>=2014 & EN4.years<=2019),3,'omitnan');

% Contour formation zones
[C,h]=m_contourf(EN4.lon_grid,EN4.lat_grid,EN4.thick_change,0:25:250,'LineColor','none');
cm=cmocean('amp');cm(1,:)=[1,1,1];%;1,1,1;1,1,1];
colormap(ax1,cm);
[ax,h]=m_contfbar(ax1,[0.3 0.80],1.05,C,h,'axfrac',.03,'FontSize',8,...
    'fontweight','bold','endpiece','no');
text(0.5,1,{'Thickness';'loss (m)'},'rotation',0,'color','k','FontSize',8,...
    'units','normalized','fontweight','bold');
% Contour mean thickness
[C2,h2]=m_contour(EN4.lon_grid,EN4.lat_grid,mean(EN4.thick_grid(:,:,EN4.years<2014),...
    3,'omitnan'));
clabel(C2,h2,'manual','FontSize',7,'Color','k');

m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','FontSize',8,'box','on',...
    'linest',':','color',[0.7,0.7,0.7]);
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','FontSize',8,'box','on');

m_line([EN4.lon_grid(1,:) EN4.lon_grid(:,end)' fliplr(EN4.lon_grid(1,:)) EN4.lon_grid(:,1)'],...
    [EN4.lat_grid(end,:) fliplr(EN4.lat_grid(:,end)') EN4.lat_grid(1,:) EN4.lat_grid(:,end)'],...
    'color','k','linestyle','-.','linewidth',1.5);
m_line(ones(1,length(23:36))*-68,23:36,'color','k','linestyle','-','linewidth',1.5);
[x,y]=meshgrid(-75:-45,30:41);
m_line([x(1,:) x(:,end)' fliplr(x(1,:)) x(:,1)'],...
    [y(end,:) fliplr(y(:,end)') y(1,:) y(:,end)'],...
    'color','k','linestyle','--','linewidth',1.5);

cm=linspecer(length(stations.coords_list(:,1))+1);
cm(1,:)=[];
for i=1:length(stations.coords_list(:,1))
    m_scatter(stations.coords_list(i,1),stations.coords_list(i,2),30,...
        'filled','markerfacecolor',cm(i,:),'markeredgecolor','k');
end
% m_scatter(-64.16666,31.666,25,'k','filled');
m_text(-64.16666+0.5,31.666-0.1,'BATS','color','k','fontweight','bold','Fontsize',7);
text(0.95,0.9,'a)','Units','normalized','FontSize',8,'FontWeight','bold');

axpos=linspace(0.875,0.03,length(stations.coords_list(:,1)));

txtlst={'b)';'c)';'d)';'e)';'f)';'g)';'h)';'i)'};
for i=1:length(axpos)-1
    ax2(i)=axes('position',[0.675 axpos(i) 0.3 0.1]);
    hold on
    ax2(i).XAxis.Visible='off';
    ax2(i).YAxis.FontSize=7;
    ax2(i).Color='none';
    text(0.05,0.1,txtlst{i},'Units','normalized','Fontsize',7,'FontWeight','bold');
end
ax2(i+1)=axes('position',[0.675 axpos(end) 0.3 0.1]);
hold on
ax2(i+1).XAxis.FontSize=7;
ax2(i+1).YAxis.FontSize=7;
ax2(i+1).Color='none';
text(0.05,0.1,txtlst{i+1},'Units','normalized','Fontsize',7,'FontWeight','bold');


for i=1:length(stations.coords_list(:,1))
    axes(ax2(i));
    tmp_thick=stations.(['Sthick' num2str(i)]);
    tmp_ORA=stations.(['ORA_' num2str(i)]);
    tmp_EN4=stations.(['EN4_' num2str(i)]);
    
    O_l(i)=plot(ax2(i),ORA.years,tmp_ORA,'color',rgb('light red'),...
        'linewidth',2);
    E_l(i)=plot(ax2(i),EN4.years,tmp_EN4,'color',rgb('light blue'),...
        'linewidth',2);
    s(i)=scatter(ax2(i),tmp_thick(:,2),tmp_thick(:,1),15,'filled',...
        'markerfacecolor',cm(i,:),'markeredgecolor','k','markerfacealpha',0.5);
    axis tight
    ylim([0 400]);
    ax2(i).YGrid='on';
    ax2(i).GridAlpha = 0.1;
end    
axes(ax2(4))
ylabel('Thickness (m)','Fontweight','bold','Fontsize',7);

tmp=1955:2018;
weak_mwt=tmp>=1971 & tmp<=1975;
weaker_mwt=tmp>=2014 & tmp<=2018;

iax=axes('position',[0.675 axpos(end) 0.3 axpos(1)+0.075]);
plot(1955:2018,NaN(length(1955:2018),1));
yl=ylim;
yl(1)=yl(1)+(max(yl)/1000);  yl(2)=yl(2)-(max(yl)/1000);
f1=patch([min(tmp(weak_mwt)) max(tmp(weak_mwt)) max(tmp(weak_mwt)) min(tmp(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp(weaker_mwt)) max(tmp(weaker_mwt)) max(tmp(weaker_mwt)) min(tmp(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
iax.XAxis.Visible='off';
iax.YAxis.Visible='off';
iax.XGrid='on';
iax.XMinorGrid='on';
iax.MinorGridLineStyle='-';
iax.GridAlpha = 0.075;
iax.MinorGridAlpha = 0.04;
% y=legend(l,{'BATS/HS','Argo','EN4g10','ORAS4'},'Fontsize',7,'Location','best');
% title(y,'Dataset');
uistack(iax,'bottom');
linkaxes([ax2,iax],'x')

load('EN4g10_data_V2.mat','dens_trans','time');
time=unique(time);
% [xn,yn]=meshgrid(0:2:600,20:1:43);

xn=repmat([0:2:600]',1,length(20:1:43));
yn=repmat(20:43,length(0:2:600),1);

dens_change=mean(dens_trans(:,:,time>=2014 & time<2019)-...
    mean(dens_trans(:,:,time<2014),3,'omitnan'),...
    3,'omitnan');

ax3=axes('position',[0.062 0.06 0.55 0.275]);
% ax3=axes('position',[0.072 0.225 0.525 0.15]);
[C,h]=contourf(yn,xn,dens_change,...
    -0.35:0.05:0.1,'linecolor','none');
colormap(ax3,cmocean('-balance','pivot',0));
[ax,h]=m_contfbar(ax3,[0.25 1],1.05,C,h,'axfrac',.03,'FontSize',8,...
        'fontweight','bold','endpiece','no');
xlim([23 36]);
xticks(24:2:36);
xticklabels(gca,{'24^\circ' '26^\circ' '28^\circ' '30^\circ' '32^\circ' '34^\circ' '36^\circ'}); 
set(gca,'ydir','reverse','fontsize',8);
xlabel('Latitude','Fontweight','bold','fontsize',8);
ylabel('Depth (m)','Fontweight','bold','fontsize',8);

ax4=axes('position',[0.062 0.06 0.55 0.275]);
contour(yn,xn,...
    mean(dens_trans(:,:,time<2014),3,'omitnan'),[26.3 26.5],'Color','k');
hold on
contour(yn,xn,...
    mean(dens_trans(:,:,time>2014 & time<2019),3,'omitnan'),...
    [26.3 26.5],'Color','k','linestyle','--');
text(-0.075,1,'j)','Units','normalized','FontSize',8,'FontWeight','bold');
text(0.5,1,{'\sigma_\theta';'change (kg m^{-3})'},'rotation',0,'color','k','FontSize',8,...
    'units','normalized','fontweight','bold');
xlim([23 36]);
set(gca,'ydir','reverse');
ax4.Visible = 'off';
ax4.XTick = [];
ax4.YTick = [];
% contourf(flipud(rot90(yn,3)),flipud(rot90(xn,3)),dens_trans(:,:,100));
linkaxes([ax3,ax4])


%%
hold on
plot(ORA.years,NaN(length(ORA.years),1));
box on
hold on
set(gca,'fontsize',7);
plot(ORA.years,smooth(ORA.outcrop_ts,3),'linewidth',2,...
    'color',rgb('light red'));
plot(EN4.years,smooth(EN4.outcrop_ts,3),'linewidth',2,...
    'color',rgb('light blue'));
yl=ylim;
yl(1)=yl(1)+(max(yl)/1000);  yl(2)=yl(2)-(max(yl)/1000);
f1=patch([min(tmp(weak_mwt)) max(tmp(weak_mwt)) max(tmp(weak_mwt)) min(tmp(weak_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
f2=patch([min(tmp(weaker_mwt)) max(tmp(weaker_mwt)) max(tmp(weaker_mwt)) min(tmp(weaker_mwt))],...
    [yl(1) yl(1) yl(2) yl(2)],rgb('very light blue'),'linestyle','none');
uistack(f1,'bottom');uistack(f2,'bottom');
yyaxis right
plot(BATS.temp_thick_annual(:,2),smooth(BATS.temp_thick_annual(:,1),3),':',...
    'linewidth',1,'color',rgb('black'));
axis tight
set(gca,'Ycolor','k');
ylabel('Thickness (m)','fontsize',8);

yyaxis left
axis tight
yl=ylim;ylim([0 max(ylim)]);
% xlim([min
grid on
ax3.GridAlpha = 0.1;
linkaxes([iax ax2 ax3],'x');
ylabel('Outcropping Volume (Svy)','fontsize',8);

%%
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure_formation.pdf -dpdf -painters
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure_formation.eps -deps 
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure_formation.png -dpng -painters
