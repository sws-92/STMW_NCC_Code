%% STMW Statistics for the paper
addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))
clear
%% BATS Data
BATS=load('20191010_extended_dataset');
% Create annual BATS values
BATS.temp_thick_annual=[];
BATS.pv_min_annual=[];
BATS.pv_props_annual=[];
BATS.intensity_annual=[];
BATS.heat_content_annual=[];
BATS.stmw_heat_content_annual=[];
BATS.years=1955:2018;

count=0;
count2=0;
for i=BATS.years(1):BATS.years(end)
        
    % isotherm thickness
    tmp_data=BATS.stevens_thickness_final(BATS.time_final>=i &...
        BATS.time_final<i+1);
    annuali=nanmean(tmp_data);
    BATS.temp_thick_annual=[BATS.temp_thick_annual;annuali i];
     
    % Intensity
    tmp_data=BATS.worthington_intensity_final(BATS.time_final>i &...
        BATS.time_final<i+1);
    intensi=nanmean(tmp_data,1);
    BATS.intensity_annual=[BATS.intensity_annual;intensi i];
   
    % PV Min. 
    tmp_data=BATS.low_pv_final((BATS.time_final>i & BATS.time_final<i+1),1);
    var_annual=nanmean(tmp_data,1);
    BATS.pv_min_annual=[BATS.pv_min_annual;var_annual i];   
    
   % PV props. 
    tmp_data=BATS.pv_data_final((BATS.time_final>i & BATS.time_final<i+1),1);
    var_annual=nanmean(tmp_data,1);
    BATS.pv_props_annual=[BATS.pv_props_annual;var_annual i];  
    
    % Heat content
    tmp_data=BATS.heat_content_final((BATS.time_final>i &...
        BATS.time_final<i+1),:);
    var_annual=nanmean(tmp_data);
    BATS.heat_content_annual=[BATS.heat_content_annual;var_annual i]; 
    
    % STMW heat content 
    tmp_data=BATS.stmw_heat_content_final((BATS.time_final>i &...
        BATS.time_final<i+1),:);
    var_annual=nanmean(tmp_data,1);
    BATS.stmw_heat_content_annual=[BATS.stmw_heat_content_annual;var_annual i]; 
end

%% Print all table values to file
HS=load('20191010_extended_dataset');
HS.pv_data_final(HS.pv_data_final(:,1)<2,1)=NaN;
fid=fopen('STMW_table.txt','w');
fprintf(fid,'%-31s%-31s%-41s%-41s%-41s%-31s%-31s%-31s%-31s\r\n',...
    'Parameter','Period','Mean (95% C.L.)','Change',...
    'Slope and std. error','n','r^2','p','%');

masks=[HS.time_final>1900 & HS.time_final<1988 HS.time_final>1988 & HS.time_final<2019 ...
    HS.time_final>2010 & HS.time_final<2019  HS.time_final>1900 & HS.time_final<2019 ...
    HS.time_final>1955 & HS.time_final<2013 HS.time_final>1967 & HS.time_final<1976];

% q=input('Do you want to refresh the combnk number? y/n','s');
% if strcmp(q,'y')
%     msksm=sum(masks,1);
%     numcomb=combnk(max(msksum));
%     save('numcomb.mat','numcomb');
% end

qs=struct('thick',HS.stevens_thickness_final,'intens',HS.worthington_intensity_final,...
    'pv',HS.low_pv_final(:,1));%,'temp',HS.pv_data_final(:,1),'surf_enth',HS.heat_content_final,'stmw_enth',HS.stmw_heat_content_final);
myfields=fieldnames(qs);

for i=1:size(masks,2)
    for ii=1:length(myfields)
        idx=logical(~isnan(qs.(myfields{ii})).*masks(:,i));
        tmp_prop=qs.(myfields{ii})(idx);
        tmp_time=HS.time_final(idx);
        
        [coeffs(1),coeffs(2)]=tsreg(tmp_time,tmp_prop,sum(idx));
%         coeffs=polyfit(tmp_time,tmp_prop,1);
%         fittedx=linspace(min(tmp_time),max(tmp_time),sum(idx));
        %         fittedy=polyval(coeffs,fittedx);
        fittedy=polyval(coeffs,tmp_time);
        r=corr(fittedy,tmp_prop)^2;
        err=sqrt((sum((fittedy-tmp_prop).^2))/(sum(idx)-2))*2/...
            sqrt(sum((tmp_time-nanmean(tmp_time)).^2));
        change(1)=(coeffs(1)-err)*length(unique(floor(tmp_time)));
        change(2)=(coeffs(1)+err)*length(unique(floor(tmp_time)));
        pct(1)=(((coeffs(1)+err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        pct(2)=(((coeffs(1)-err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        [~,sig]=Mann_Kendall(tmp_prop,0.05);
        SEM=nanstd(tmp_prop)/sqrt(length(tmp_prop)); % Standard Error
        ts=tinv([0.025  0.975],length(tmp_prop)-1); % T-Score
        cl=nanmean(tmp_prop)+ts*SEM; % Confidence Intervals

        %         if i==3 && ii==1
        %         keyboard
%         end
        fprintf(fid,'%-30s %-4.0f-%-25.0f %-10.4e (%-10.4e-%-16.4e) %-15.4e-%-25.4e %-15.4e+-%-25.4e %-30.0f %-30.2f %-30.2f %-4.0f-%-25.0f\r\n',...
            myfields{ii},floor(min(tmp_time)),ceil(max(tmp_time)),mean(tmp_prop),...
            cl(1),cl(2),change(1),change(2),coeffs(1),err,...
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

masks=[HS.time_final>1900 & HS.time_final<1988 HS.time_final>1988 & HS.time_final<2019 ...
    HS.time_final>2010 & HS.time_final<2019  HS.time_final>1900 & HS.time_final<2019 ...
    HS.time_final>1955 & HS.time_final<2013 HS.time_final>1966 & HS.time_final<1974];

qs=struct('temp',HS.pv_data_final(:,1),'STMWenth',HS.stmw_heat_content_final,...
    'Surfenth',HS.heat_content_final);%,'temp',HS.pv_data_final(:,1),'surf_enth',HS.heat_content_final,'stmw_enth',HS.stmw_heat_content_final);
myfields=fieldnames(qs);

for i=1:size(masks,2)
    for ii=1:length(myfields)
        idx=logical(~isnan(qs.(myfields{ii})).*masks(:,i));
        tmp_prop=inpaint_nans(qs.(myfields{ii})(idx));
        tmp_time=HS.time_final(idx);
        
        [coeffs(1),coeffs(2)]=tsreg(tmp_time,tmp_prop,sum(idx));
%         coeffs=polyfit(tmp_time,tmp_prop,1);
%         fittedx=linspace(min(tmp_time),max(tmp_time),sum(idx));
        %         fittedy=polyval(coeffs,fittedx);
        fittedy=polyval(coeffs,tmp_time);
        r=corr(fittedy,tmp_prop)^2;
        err=sqrt((sum((fittedy-tmp_prop).^2))/(sum(idx)-2))*2/...
            sqrt(sum((tmp_time-nanmean(tmp_time)).^2));
        change(1)=(coeffs(1)-err)*length(unique(floor(tmp_time)));
        change(2)=(coeffs(1)+err)*length(unique(floor(tmp_time)));
        pct(1)=(((coeffs(1)+err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        pct(2)=(((coeffs(1)-err)*length(unique(floor(tmp_time))))/fittedy(1))*100;
        [~,sig]=Mann_Kendall(tmp_prop,0.05);
        SEM=nanstd(tmp_prop)/sqrt(length(tmp_prop)); % Standard Error
        ts=tinv([0.025  0.975],length(tmp_prop)-1); % T-Score
        cl=nanmean(tmp_prop)+ts*SEM; % Confidence Intervals
        
        fprintf(fid,'%-30s %-4.0f-%-25.0f %-10.4e+-%-30.4e %-15.4e-%-25.4e %-15.4e+-%-25.4e %-30.0f %-30.2f %-30.2f %-4.0f-%-25.0f\r\n',...
            myfields{ii},floor(min(tmp_time)),ceil(max(tmp_time)),mean(tmp_prop),...
            std(tmp_prop),change(1),change(2),coeffs(1),err,...
            sum(idx),r,sig,pct(1),pct(2));
        
        if strcmp(myfields{ii},'temp') && sum(idx)>2000 && sum(idx)<2500
            fprintf(['2010-2019 per decade temp. change:\n %3.2f - %3.2f \n'],...
                (coeffs(1)-err)*10,(coeffs(1)+err)*10);
        elseif strcmp(myfields{ii},'temp') && sum(idx)>5000
            fprintf(['Mean STMW temp: %3.2f \n'],...
                nanmean(tmp_prop));
            fprintf(['Mean 2017 STMW temp: %3.2f \n'],...
                nanmean(tmp_prop(tmp_time>2017 & tmp_time<2018)));
            fprintf(['Mean 2018 STMW temp: %3.2f \n'],...
                nanmean(tmp_prop(tmp_time>2018 & tmp_time<2019)));
        end
    end
    
    fprintf(fid,'\r\n');
end
fclose(fid);

recent_mwt=HS.time_final>2014 & HS.time_final<2019;
pct=((sum(HS.stevens_thickness_final(recent_mwt)==0))/...
    length(HS.stevens_thickness_final(recent_mwt)))*100;
disp([num2str(pct),'% of casts =0m thickness (2014-2019)']);

tmp_prop=inpaint_nans(BATS.pv_props_annual(:,1));
tmp_time=[1955:2018]';
[coeffs(1),coeffs(2)]=tsreg(tmp_time,tmp_prop,...
    length(tmp_time));
fittedy=polyval(coeffs,tmp_time);
err=sqrt((sum((fittedy-tmp_prop).^2))/(length(tmp_time)-2))*1.96/...
    sqrt(sum((tmp_time-nanmean(tmp_time)).^2));
fprintf(['Overall per century temp. change:\n %3.2f - %3.2f \n'],...
    (coeffs(1)-err)*100,(coeffs(1)+err)*100);

%% Formation Stats
ORA=load('ORAS4_all.mat');
tmp=load('ORAS4_data.mat','time');ORA.all_time=tmp.time;
EN4=load('EN4g10_all.mat');

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
count=0;

for i=1:length(EN4.years)
    idx=EN4.time>EN4.years(i) & EN4.time<EN4.years(i)+1;
    EN4.annual_thick(:,:,i)=mean(EN4.thick_grid(:,:,idx),3,'omitnan');
    EN4.annual_intens(:,:,i)=mean(EN4.intensity_grid(:,:,idx),3,'omitnan');
    EN4.annual_pv(:,:,i)=mean(EN4.PV_grid(:,:,idx),3,'omitnan');
    EN4.annual_stemp(:,:,i)=mean(EN4.stemp_grid(:,:,idx),3,'omitnan');
    EN4.outcrop_ts(i)=sum(sum(EN4.outcrop_grid(:,:,i),1),2)/(1e6*60*60*24*365);
    EN4.hc_ts(i)=sum(sum(sum(EN4.hc_grid(:,:,idx),1,'omitnan'),2,'omitnan'));

    count=count+1;
    idx=floor(EN4.time)==EN4.years(i);
    EN4.stemp_yr(:,:,count)=mean(EN4.stemp_grid(:,:,idx),3,'omitnan');
    
    % interplate lat to 18 degree temp. surface
    for j=1:length(EN4.lon_grid(1,:))
        if EN4.lon_grid(1,j)>-75 && EN4.lon_grid(1,j)<-45 && sum(idx)~=0
            tmp=EN4.stemp_yr(~isnan(EN4.stemp_yr(:,j,count)),j,count);
            [tmp,uniq_idx]=unique(tmp);
            EN4.OC_lat(j,count)=interp1(tmp,EN4.lat_grid(uniq_idx,1),18,'linear');
        else
            EN4.OC_lat(j,count)=NaN;
        end
    end   
    EN4.OC_position(count)=nanmean(EN4.OC_lat(:,count));

end


count=0;
for i=1:length(ORA.years)
    idx=ORA.time>ORA.years(i) & ORA.time<ORA.years(i)+1;
    ORA.annual_thick(:,:,i)=mean(ORA.thick_grid(:,:,idx),3,'omitnan');
    ORA.annual_intens(:,:,i)=mean(ORA.intensity_grid(:,:,idx),3,'omitnan');
    ORA.annual_pv(:,:,i)=mean(ORA.PV_grid(:,:,idx),3,'omitnan');
    ORA.annual_stemp(:,:,i)=mean(ORA.stemp_grid(:,:,idx),3,'omitnan');
    ORA.outcrop_ts(i)=sum(sum(ORA.outcrop_grid(:,:,i),1),2)/(1e6*60*60*24*365);
    ORA.hc_ts(i)=sum(sum(sum(ORA.hc_grid(:,:,idx),1,'omitnan'),2,'omitnan'));
    
    count=count+1;
    idx=floor(ORA.time)==ORA.years(i);
    ORA.stemp_yr(:,:,count)=mean(ORA.stemp_grid(:,:,idx),3,'omitnan');
    
    % interplate lat to 18 degree temp. surface
    for j=1:length(ORA.lon_grid(1,:))
        if ORA.lon_grid(1,j)>-75 && ORA.lon_grid(1,j)<-45 && sum(idx)>1           
            tmp=ORA.stemp_yr(~isnan(ORA.stemp_yr(:,j,count)),j,count);
           [tmp,uniq_idx]=unique(tmp);
            ORA.OC_lat(j,count)=interp1(tmp,ORA.lat_grid(uniq_idx,1),18,'linear');
        else
            ORA.OC_lat(j,count)=NaN;
        end
    end
    ORA.OC_position(count)=nanmean(ORA.OC_lat(:,count));
end

clc
[RHO, PVAL]=corr(inpaint_nans(BATS.intensity_annual(:,1)),EN4.outcrop_ts(1:end-1)');
disp(['BATS intensity- EN4 OC volume r=',num2str(RHO)]);
[RHO, PVAL]=corr(inpaint_nans(BATS.intensity_annual(4:end-1,1)),ORA.outcrop_ts(1:end-1)');
disp(['BATS intensity- ORA OC volume r=',num2str(RHO)]);

[RHO, PVAL]=corr(inpaint_nans(BATS.temp_thick_annual(:,1)),EN4.outcrop_ts(1:end-1)');
disp(['BATS thickness- EN4 OC volume r=',num2str(RHO)]);
[RHO, PVAL]=corr(inpaint_nans(BATS.temp_thick_annual(4:end,1)),ORA.outcrop_ts(1:end)');
disp(['BATS thickness- ORA OC volume r=',num2str(RHO)]);

[r,pval]=corr(inpaint_nans(BATS.heat_content_annual(:,1)),EN4.outcrop_ts(1:end-1)');
disp(['BATS surface heat-EN4 OC volume r=',num2str(r)]);
[r,pval]=corr(inpaint_nans(BATS.heat_content_annual(4:end,1)),ORA.outcrop_ts(1:end)');
disp(['BATS surface heat-ORA OC volume r=',num2str(r)]);

[r,pval]=corr(EN4.hc_ts(1:end-1)',EN4.outcrop_ts(1:end-1)');
disp(['EN4 surface heat-EN4 OC volume r=',num2str(r)]);
[r,pval]=corr(ORA.hc_ts(1:end-1)',ORA.outcrop_ts(1:end-1)');
disp(['ORA surface heat-ORA OC volume r=',num2str(r)]);

[r,pval]=corr(inpaint_nans(BATS.intensity_annual(:,1)),EN4.OC_position(1:end-1)');
disp(['BATS intensity-EN4 OC position r=',num2str(r)]);
[r,pval]=corr(inpaint_nans(BATS.intensity_annual(4:end-1,1)),ORA.OC_position(1:end-1)');
disp(['BATS intensity-ORA OC position r=',num2str(r)]);

[r,pval]=corr(inpaint_nans(BATS.temp_thick_annual(:,1)),EN4.OC_position(1:end-1)');
disp(['BATS thickness-EN4 OC position r=',num2str(r)]);
[r,pval]=corr(inpaint_nans(BATS.temp_thick_annual(4:end-1,1)),ORA.OC_position(1:end-1)');
disp(['BATS thickness-ORA OC position r=',num2str(r)]);

[r,pval]=corr(EN4.OC_position(1:end-1)',EN4.outcrop_ts(1:end-1)');
disp(['EN4 OC position-EN4 OC volume r=',num2str(r)]);
[r,pval]=corr(ORA.OC_position(1:end-1)',ORA.outcrop_ts(1:end-1)');
disp(['ORA OC position-ORA OC volume r=',num2str(r)]);


% disp(['BATS intensity-OC position r=',num2str(r)]);
% [r,pval]=corr(OC_volume,OC_position);
% disp(['OC volume-OC position r=',num2str(r)]);

%% OC volume change
idx=EN4.years>=2010 & EN4.years<2019;
[RHO, PVAL]=Mann_Kendall(EN4.outcrop_ts(idx),0.05);
[coeffs(1),coeffs(2)]=tsreg(EN4.years(idx),EN4.outcrop_ts(idx),sum(idx));
fittedx=polyval(coeffs,EN4.years(idx));
err=sqrt((sum((fittedx-EN4.outcrop_ts(idx)).^2))/(length(EN4.years)-2))*2/...
    sqrt(sum((EN4.years(idx)-nanmean(EN4.years(idx))).^2));
change(1)=(((coeffs(1)-err)*length(EN4.years(idx)))/fittedx(1))*-100;
change(2)=(((coeffs(1)+err)*length(EN4.years(idx)))/fittedx(1))*-100;

tmpa=sprintf('%2.2f',coeffs(1));
tmpb=sprintf('%2.1f',change(1));
tmpc=sprintf('%2.1f',change(2));
tmpd=sprintf('%2.2f',err);
disp([tmpa,' +- ',tmpd,' Svy year reduction in EN4 OC volume']);
disp([tmpc,'-',tmpb,'% reduction in EN4 OC volume']);

tmpa=sprintf('%2.1f',(length(EN4.years(idx))*(-1*coeffs(1)+err)));
tmpb=sprintf('%2.1f',(length(EN4.years(idx))*(-1*coeffs(1)-err)));
disp([tmpb,' - ',tmpa,' Svy overall reduction in EN4 OC volume']);

idx=ORA.years>2009 & ORA.years<2018;
[RHO, PVAL]=Mann_Kendall(ORA.outcrop_ts(idx),0.05);
[coeffs(1),coeffs(2)]=tsreg(ORA.years(idx),ORA.outcrop_ts(idx),sum(idx));
fittedx=polyval(coeffs,ORA.years(idx));
err=sqrt((sum((fittedx-ORA.outcrop_ts(idx)).^2))/(length(ORA.years)-2))*2/...
    sqrt(sum((ORA.years(idx)-nanmean(ORA.years(idx))).^2));
change(1)=(((coeffs(1)-err)*length(ORA.years(idx)))/fittedx(1))*-100;
change(2)=(((coeffs(1)+err)*length(ORA.years(idx)))/fittedx(1))*-100;

tmpa=sprintf('%2.2f',coeffs(1));
tmpb=sprintf('%2.1f',change(1));
tmpc=sprintf('%2.1f',change(2));
tmpd=sprintf('%2.2f',err);
disp([tmpa,' +- ',tmpd,' Svy year reduction in ORA OC volume']);
disp([tmpc,'-',tmpb,'% reduction in ORA OC volume']);

tmpa=sprintf('%2.1f',(length(ORA.years(idx))*(-1*coeffs(1)+err)));
tmpb=sprintf('%2.1f',(length(ORA.years(idx))*(-1*coeffs(1)-err)));
disp([tmpb,' - ',tmpa,' Svy overall reduction in ORA OC volume']);

% Check pentads of OC position
count=0;
OC_pos_pent=NaN(length(EN4.years-5),1);
for i=EN4.years(1):EN4.years(end)-5
    count=count+1;
    time_idx=EN4.years>i & EN4.years<i+5;
    OC_pos_pent(count)=nanmean(EN4.OC_position(time_idx));
end
[tmp, idx]=max(OC_pos_pent);
disp('The pentad with the northernmost 18^oC outcropping was:');
disp([num2str(EN4.years(idx)),'-',num2str(EN4.years(idx)+4),' at ',num2str(tmp),'^oN']);

% [RHO, PVAL]=corr(BATS.intensity_annual(5:end,1),OC_position(1:end-1));
% [RHO, PVAL]=corr(sat_intensity,OC_position);
% 
[RHO, PVAL]=Mann_Kendall(EN4.OC_position,0.05);


%%
% tmpa=sprintf('%2.1f',(length(years)*(-1*coeffs(1)+err))/3.15e13);
% tmpb=sprintf('%2.1f',(length(years)*(-1*coeffs(1)-err))/3.15e13);
% disp([tmpb,' - ',tmpa,' Svy reduction in OC volume']);

[coeffs(1),coeffs(2)]=tsreg(years,sat_intensity,1);
fittedx=polyval(coeffs,years);
[RHO, PVAL]=corr(BATS.intensity_annual(5:end,1),sat_intensity(1:end-1));
[RHO, PVAL]=Mann_Kendall(sat_intensity,0.05);
% disp([num2str(((max(fittedx)-min(fittedx))/max(fittedx))*100),'% reduction in OC intensity']);

% coeffs=polyfit(years,formation_pv_min_temp(:,1),1);


