%% PLOTS for paper

addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))

lnes=lines;

% if ~exist('stevens_thickness_final','var')
load('20191010_extended_dataset.mat');
% end
%% Window averaging

% pv_cont_final=pv_cont(:,(time>0));
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
timecount=1955:1/2:max(time_final);
depth=5:5:600;
tmonth=zeros(length(depth),length(timecount));
pmonth=zeros(length(depth),length(timecount));
tgmonth=zeros(length(depth),length(timecount));
wlimmonth=zeros(length(timecount),2);
rmonth=zeros(length(depth),length(timecount));
rgmonth=zeros(length(depth),length(timecount));

for i=1:length(timecount)
    ctime=timecount(i);
    tmonth(:,i)=nanmean(temp_cont_final(:,(time_final>=ctime & time_final<ctime+1/2)),2);
    pmonth(:,i)=nanmean(pv_cont_final(:,(time_final>=ctime & time_final<ctime+1/2)),2);
    tgmonth(:,i)=nanmean(tgrad_cont_final(:,(time_final>=ctime & time_final<ctime+1/2)),2)*50;
    wlimmonth(i,1)=nanmean(stevens_limits_final((time_final>=ctime & time_final<ctime+1/2),1));
    wlimmonth(i,2)=nanmean(stevens_limits_final((time_final>=ctime & time_final<ctime+1/2),2));
    rmonth(:,i)=nanmean(dens_cont_final(:,(time_final>=ctime & time_final<ctime+1/2)),2);
    rgmonth(:,i)=nanmean(d_grad_cont_final(:,(time_final>=ctime & time_final<ctime+1/2)),2)*50;
end

% 
% timecount(isnan(wlimmonth(:,1)))=[];
% tmonth(:,(isnan(tmonth(100,:))))=[];
% tmonth=inpaint_nans(tmonth);
% pmonth(:,(isnan(pmonth(100,:))))=[];
% pmonth=inpaint_nans(pmonth);
% wlimmonth((isnan(wlimmonth(:,1))),:)=[];
% elimmonth((isnan(elimmonth(:,1))),:)=[];

%% Results plots
% Isotherm
f1=figure('units','centimeters','outerposition',[0.01 0.01 8.9 12],'color','w');
ax1=axes('position',[.1 .7 .70 .25]);
[C,h]=contourf(timecount,depth,tmonth,'LineColor','none','LevelList',[13:0.5:17 18 19:.5:23]);
colormap(ax1,cmocean('thermal','pivot',19));
a=m_contfbar(ax1,1.05,[0.0 1],C,13:23,'axfrac',.035,'fontsize',7,'endpiece','no');
% caxis([14 21])
ylabel(a,'\Theta (^oC)','FontSize',7,'rotation',0,'FontWeight','bold');
set(gca,'ydir','reverse','FontSize',7);
hold on
[C,h]=contour(timecount,depth,tmonth,[17 19],'color',...
    [0.9 0.9 0.9],'linewidth',1);
[C,h]=contour(timecount,depth,tmonth,[18 18],'color',...
    [0.9 0.9 0.9],'linewidth',1,'linestyle',':');
ylim([50 550]);
text(-0.1,1,'a)','Units','normalized','FontSize',7,'FontWeight','bold')

% PV
ax2=axes('position',[.1 .4 .70 .25]);
[C,h]=contourf(timecount,depth,pmonth,'LineColor','none','LevelList',...
    [1e-11:5e-12:2e-10]);
colormap(ax2,m_colmap('diverging'));
a=m_contfbar(ax2,1.05,[0.0 1],C,[1e-11:20e-12:2e-10],'axfrac',.035,'fontsize',7,'endpiece','no',...
    'levels','set');
ylabel('Depth (m)','FontSize',7,'FontWeight','bold');
ylabel(a,{'PV';'(s^{-1} m^{-1})'},'FontSize',7,'rotation',0,'FontWeight','bold');
set(gca,'ydir','reverse','FontSize',7);
hold on
[C,h]=contour(timecount,depth,pmonth,[1e-10 1e-10],'color',...
    [0.4 0.4 0.4],'linewidth',1);
% plot(time_final,movmean(low_pv_final(:,2),20),'k','linewidth',2);
ylim([50 550]);
text(-0.1,1.1,'b)','Units','normalized','FontSize',7,'FontWeight','bold')

% Isopycnals
ax3=axes('position',[.1 .1 .70 .25]);
[C,h]=contourf(timecount,depth,rmonth,[24:.2:26.2 26.33 26.42 26.51 26.6:.2:27],...
    'LineColor','none');
colormap(ax3,cmocean('curl','pivot',26.33));
a=m_contfbar(ax3,1.05,[0.0 1],C,24:.4:27,'axfrac',.035,'fontsize',7,'endpiece','no',...
    'levels','set');
ylabel(a,{'\sigma_{\theta}';'(kg m^{-3})'},'FontSize',7,'rotation',0,'FontWeight','bold');
set(gca,'ydir','reverse','FontSize',7);
hold on
[C,h]=contour(timecount,depth,rmonth,[26.33 26.51],'color',...
    [0.4 0.4 0.4],'linewidth',1);
[C,h]=contour(timecount,depth,rmonth,[26.42 26.42],'color',...
    [0.4 0.4 0.4],'linewidth',1,'linestyle',':');
% plot(time_final,movmean(low_pv_final(:,2),20),'k','linewidth',2);
ylim([50 550]);
text(-0.1,1.1,'c)','Units','normalized','FontSize',7,'FontWeight','bold')

%%
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure_contour.pdf
export_fig C:\Users\samst\Dropbox\UBC\STMW_project\STMW_NCC_Review1\Figures\Figure_contour.eps -deps
