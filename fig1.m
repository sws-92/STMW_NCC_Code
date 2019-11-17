%% Create Fig 1
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

%% Fig 1
all_temp=mean(EN4.wtemp_grid(:,:,EN4.years<2013),3,'omitnan');
rec_temp=mean(EN4.wtemp_grid(:,:,EN4.years>=2014 & EN4.years<=2018),3,'omitnan');

% Set up figure
f1=figure('units','centimeters','outerposition',[0.01 0.01 18 13],'color','w');
ax1=axes('position',[-0.15 0.2 1 0.65]);
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
[ax,h]=m_contfbar([0.22 0.7],1.05,C,h,'axfrac',.03,'FontSize',8,'fontweight','bold');
ylabel(ax,{'SST';'(^oC)'},'rotation',0,'color','k');
[C,c]=m_contour(EN4.lon_grid,EN4.lat_grid,all_temp,[17 19],'w--','linewidth',1.5);
[B,b]=m_contour(EN4.lon_grid,EN4.lat_grid,all_temp,[18 18],'w-.','linewidth',1.5);
[F,f]=m_contour(EN4.lon_grid,EN4.lat_grid,rec_temp,[17 19],'k--','linewidth',1);
[F,f]=m_contour(EN4.lon_grid,EN4.lat_grid,rec_temp,[18 18],'k-.','linewidth',1);

m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','FontSize',8,'box','on');
m_line([EN4.lon_grid(1,:) EN4.lon_grid(:,end)' fliplr(EN4.lon_grid(1,:)) EN4.lon_grid(:,1)'],...
    [EN4.lat_grid(end,:) fliplr(EN4.lat_grid(:,end)') EN4.lat_grid(1,:) EN4.lat_grid(:,end)'],...
    'color','k','linestyle','-.','linewidth',1.5);
[x,y]=meshgrid(-75:-45,30:41);
m_line([x(1,:) x(:,end)' fliplr(x(1,:)) x(:,1)'],...
    [y(end,:) fliplr(y(:,end)') y(1,:) y(:,end)'],...
    'color','k','linestyle','--','linewidth',1.5);
m_scatter(-64.16666,31.666,25,'k','filled');
m_text(-64.16666+0.5,31.666-0.1,'BATS','color','k','fontweight','bold','FontSize',6);
m_scatter(-64.5,32.16666,25,'k','filled');
m_text(-64.5+0.5,32.16666+0.1,'Hydrostation ''S''','color','k','fontweight','bold','FontSize',6);
% text(0.9,0.9,'a)','Units','normalized','FontSize',8,'FontWeight','bold');

m_scatter(stations.coords_list(:,1),stations.coords_list(:,2),15,...
    'k','filled');
m_text(stations.coords_list([1 3:8],1)+0.3,stations.coords_list([1 3:8],2),...
    {'BV1','BV3','BV4','BV5','BV6','BV7','BV8'},...
    'color','k','FontSize',6);
m_text(-38,42,'a)','FontSize',7,'fontweight','bold');

%% FIGURE- PV profiles
if ~exist('HS','var')
    HS=load('20191010_extended_dataset');
end
timeyear=1955:2017;
[hc,edges]=histcounts(floor(HS.time_final),1955:2017);

% f1=figure('units','centimeters','outerposition',[0.01 0.01 6 13]);
% set(gcf,'color','w');
% subplot(3,1,1)

axes('position',[0.7 0.81 0.25 0.15]);
histogram(floor(HS.time_final),1955:2018,'facecolor',rgb('grey'));
set(gca,'yscale','log','ytick',[1 10 50 100 500]);
hold on
% Have to add these in as matlab doesnt plot them...
patch([1974 1975 1975 1974],[0 0 hc(edges==1974) hc(edges==1974)],...
    rgb('grey'),'edgecolor','k','linewidth',1);
patch([1981 1982 1982 1981],[0 0 hc(edges==1981) hc(edges==1981)],...
    rgb('grey'),'edgecolor','k','linewidth',1);
axis tight
ylim([0 600]);
ylabel({'Number of';'profiles'},'FontSize',7,'FontWeight','bold');
text(0.02,0.9,'b)','Units','normalized','FontSize',7,'FontWeight','bold')
% xlabel('Year','FontSize',8,'FontWeight','bold');
xtick([1960:15:2020]);
grid on
set(gca,'GridLineStyle','--','tickdir','out','fontsize',7);

mn_profile=nanmean(HS.pv_cont_final,2);
pent_profile=nanmean(HS.pv_cont_final(:,HS.time_final>=2014 &...
    HS.time_final<=2018),2);
dpth_prf=5:5:600;
idxf=find(mn_profile<1e-10,1,'first');
idxl=find(mn_profile<1e-10,1,'last');

axes('position',[0.7 0.2 0.25 0.55]);
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
set(gca,'ydir','reverse','fontsize',7,'tickdir','out');
axis tight
ylim([100 600]);
xlim([0 5e-10]);
text(0.9,0.05,'c)','Units','normalized','FontSize',7,'FontWeight','bold')
ylabel('Depth (m)','FontSize',7,'FontWeight','bold');
xlabel('PV (m^{-1} s^{-1})','FontSize',7,'FontWeight','bold');
uistack(p1,'bottom');
uistack(p2,'bottom');uistack(p2,'up');
uistack(p3,'bottom');uistack(p3,'up',2);
% text(3.5e-10,275,'1955-2017','color','b','fontweight','bold','fontsize',8);
annotation('textbox',[0.85 0.57 0.01 0.01],'String','1955-2018','Fitboxtotext','on',...
    'color','b','fontweight','bold','FontSize',7,'backgroundcolor','w');
annotation('textbox',[0.85 0.46 0.01 0.01],'String','2014-2018','Fitboxtotext','on',...
    'color','r','fontweight','bold','FontSize',7,'backgroundcolor','w');
box on
%%
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\map.eps -deps

