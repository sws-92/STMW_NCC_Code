%% Read EN4g10 data
clear
addpath(genpath('/ocean/sstevens/IW_project'))
 
dir_struct1=dir('/ocean/sstevens/EN4g10/data/*.nc');
 
lat=ncread('EN.4.2.1.f.analysis.g10.196001.nc','lat');
lon=ncread('EN.4.2.1.f.analysis.g10.196001.nc','lon');
xcut=-80:-40;ycut=20:1:43;
% [xgrid,ygrid]=meshgrid(xcut,ycut);
% xgrid=rot90(xgrid);ygrid=rot90(ygrid);dgrid=double(rot90(dgrid));
depth=ncread('EN.4.2.1.f.analysis.g10.196001.nc','depth');
 
lat_idx=find(lat>=20 & lat<=43);
lon_idx=find(lon>=100+180 & lon<=140+180);
depth_idx=depth<700;
years=1955:2018;
temp_grid=zeros(length(0:2:600),length(years)*12*...
    length(lat_idx)*length(lon_idx));
BATS_temp_month=[];
BATS_salt_month=[];
salt_grid=zeros(length(0:2:600),length(years)*12*...
    length(lat_idx)*length(lon_idx));
lat_vec=zeros(1,length(years)*12*length(lat_idx)*length(lon_idx));
lon_vec=zeros(1,length(years)*12*length(lat_idx)*length(lon_idx));
time=zeros(1,length(years)*12*length(lat_idx)*length(lon_idx));
[~,b_lon_idx]=min(abs(xcut--64.166));
[~,b_lat_idx]=min(abs(ycut-31.666));
 
count=0;
count2=0;
count3=0;
m_count=0;
  
% m_proj('equidistant','lon',[-80 -15],'lat',[25 50]);
% % m_proj('equidistant','lon',[-180 180],'lat',[-90 90]);
% m_coast('patch',[0.7 0.7 0.7]);
% hold on
% % m_contourf(x,y,squeeze(salt(:,:,1,1)));
% % m_scatter(x_cut(:),y_cut(:));
% % m_contourf(xcut,ycut,rot90(salt(lon_idx,lat_idx,1)));
% m_contourf(xcut,fliplr(ycut),rot90(squeeze(temp(lon_idx,lat_idx,1))));

for i=1:length(dir_struct1)
    
    ds=dir_struct1(i);
    disp(ds.name)
    salt=ncread(ds.name,'salinity');
    salt_cut=salt(lon_idx,lat_idx,:);
    tmp_t=datestr(double(ncread(ds.name,'time'))+datenum(1800,01,01),23);
    
    tmp_tstr=tmp_t;
    tmp_tstr2(1,1)=str2num([tmp_tstr(7:10) tmp_tstr(1:2) tmp_tstr(4:5)]);
    tmp_tstr2(1,2)=0000;
    tmp_t_dec=cnv_date(tmp_tstr2);
    
    BATS_salt_month=[BATS_salt_month interp1(depth,...
        squeeze(salt_cut(b_lon_idx,b_lat_idx,:)),[0:2:600]')];
    
    for ii=1:length(lon_idx)
        for iii=1:length(lat_idx)
            count2=count2+1;
            salt_grid(:,count2)=interp1(depth(depth_idx),...
                squeeze(salt_cut(ii,iii,depth_idx)),[0:2:600]');
            lon_vec(count2)=xcut(ii);
            lat_vec(count2)=ycut(iii);
            time(count2)=tmp_t_dec;
        end
    end
    
    temp=ncread(ds.name,'temperature');
    temp_cut=temp(lon_idx,lat_idx,:);
    
    BATS_temp_month=[BATS_temp_month interp1(depth,...
        squeeze(temp_cut(b_lon_idx,b_lat_idx,:)),[0:2:600]')];
    
    for ii=1:length(lon_idx)
        for iii=1:length(lat_idx)
            count3=count3+1;
            temp_grid(:,count3)=interp1(depth(depth_idx),...
                squeeze(temp_cut(ii,iii,depth_idx)),[0:2:600]');
        end
    end  
end
 
idx=salt_grid(1,:)==0;time(idx)=[];lat_vec(idx)=[];lon_vec(idx)=[];
salt_grid(:,idx)=[];temp_grid(:,idx)=[];
month=unique(time); BATS_salt_month=double(BATS_salt_month);
BATS_temp_month=double(BATS_temp_month);
 
clearvars -except lat_vec lon_vec time temp_grid salt_grid...
    BATS_temp_month BATS_salt_month
save('EN4g10_data_V2');

%% EN4g10- BATS
clear
addpath(genpath('C:\Users\samst\Documents'))
load EN4g10_data.mat

worthington_intensity=zeros(size(BATS_temp_month,2),1);
pv_data=zeros(size(BATS_temp_month,2),5);
w_grad=zeros(size(BATS_temp_month,2),1);
mld=zeros(size(BATS_temp_month,2),3);
temp_cont=zeros(600/5,size(BATS_temp_month,2));
sal_cont=zeros(600/5,size(BATS_temp_month,2));
dens_cont=zeros(600/5,size(BATS_temp_month,2));
pv_cont=zeros(600/5,size(BATS_temp_month,2));
tgrad_cont=zeros(600/5,size(BATS_temp_month,2));
heat_content=zeros(size(BATS_temp_month,2),1);
stmw_heat_content=zeros(size(BATS_temp_month,2),1);
stmw_temp=zeros(size(BATS_temp_month,2),1);
stevens_thickness=zeros(size(BATS_temp_month,2),1);
stevens_limits=zeros(size(BATS_temp_month,2),2);
outcropping=zeros(size(BATS_temp_month,2),1);
core_aou=zeros(size(BATS_temp_month,2),1);
no_layer=0;
error=[];
depth=[0:2:600]';

% Find how many profiles have EDW outcropping and run tests to ensure they
% meet EDW criteria

for i=1:size(BATS_temp_month,2)

    display(i);
    
    pres=gsw_p_from_z(depth*-1,31.5);
    
    % Calculate TEOS-10 variables
    SA=gsw_SA_from_SP(BATS_salt_month(:,i),pres,31.5,-64.5);
    CT=gsw_CT_from_pt(SA,BATS_temp_month(:,i));
    
    % Calculate density from variables
    p_dens=gsw_sigma0(SA,CT);
    p_dens=inpaint_nans(p_dens);
    
    % Calculate Worthington thickness (thickness between 17 and 19
    % degree isotherms)
    [~,idy]=min(abs(17-CT));
    mw_low=depth(idy);
    
    worthington_limits(i,2)=mw_low;
    
    [~,idx]=min(abs(19-CT));
    mw_high=depth(idx);
    worthington_limits(i,1)=mw_high;
    
    % Calculate "Stevens thickness" (thickness between 17 and 19
    % degree isotherms and PV < 1x10^-10)
    [~,q,p_ave]=sw_bfrq(BATS_salt_month(:,i),gsw_t_from_CT(SA,CT,pres),pres,31.5);
    q=interp1(p_ave,q,pres);
    
    % Filter PV signal with a moving mean to remove the high frequency variability
    q=movmean(q,35);
    msk=zeros(length(CT),1);
    msk(idx:idy)=1;
    s_msk=logical(msk.*(q<1e-10));
    
    if sum(s_msk)>2
        stevens_thickness(i)=max(depth(s_msk))-...
            min(depth(s_msk));
        stevens_limits(i,:)=[min(depth(s_msk)) ...
            max(depth(s_msk))];
    end
    
    % Calulate MLD
    try
        mld(i,1)=ra_mld(SA,CT,depth,0.5); % 1° MLD temp as per Hanawa and Hoshino (1988) and Sufa and Hanawa (1990)
        mld(i,2)=ra_mld(SA,CT,depth,1);
    catch
        mld(i,:)=NaN;
    end
    
    % Create temperature contour
    try
        temp_cont(:,i)=interp1(depth,CT,5:5:600);
    catch
        temp_cont(:,i)=NaN;
    end
    
    % Create salinity contour
    try
        sal_cont(:,i)=interp1(depth,SA,5:5:600);
    catch
        sal_cont(:,i)=NaN;
    end
    
    % Create temperature and density gradients
    CT_grad=diff(CT);
    p_dens_grad=diff(p_dens);
    
    % Create temp. gradient contour
    try
        tgrad_cont(:,i)=interp1(p_ave,CT_grad,5:5:600);
    catch
        tgrad_cont(:,i)=NaN;
    end
    
    % Calculate enthalpy
    idxh=~isnan(CT) & ~isnan(SA) &...
        ~isnan(pres) & pres<150 &...
        pres>10;
    
    % Create surface integrated heat content <150m
    heat_content(i)=trapz(pres(idxh),gsw_enthalpy(SA(idxh),...
        CT(idxh),pres(idxh)));
    
    % Create STMW integrated heat content
    if sum(s_msk)>2
        stmw_heat_content(i)=trapz(pres(s_msk),gsw_enthalpy(SA(s_msk),...
            CT(s_msk),pres(s_msk)));
        
        stmw_temp(i)=nanmean(CT(s_msk));
    else
        stmw_heat_content(i)=NaN;
        
        stmw_temp(i)=NaN;
    end
    
    % Outcropping?
    OC_idx=pres<10;
    if sum(OC_idx)>2 && nanmean(CT(OC_idx))<19
        outcropping(i)=1;
    end
    
    % Calulate core pv between 16oC isotherms
    if sum(s_msk)>2
        try
            [~,idyT]=min(abs(16-CT));[~,idxT]=min(abs(20-CT));
            low_pv(i,1)=min(q(idxT:idyT)); low_pv(i,2)=depth(q==min(q(idxT:idyT)));
            
            if low_pv(i,1)<0
                low_pv(i,1)=min(q(idx:idy)); low_pv(i,2)=depth(q==min(q(idx:idy)));
                % Extract all properties at the core PV depth level
                pv_data(i,1)=CT(q==min(q(idx:idy)));
                pv_data(i,2)=SA(q==min(q(idx:idy)));
                pv_data(i,3)=p_dens(q==min(q(idx:idy)));
                pv_data(i,4)=CT_grad(q==min(q(idx:idy)));
                pv_data(i,5)=p_dens_grad(q==min(q(idx:idy)));
                %                 core_aou(i)=tmp_aou(q==min(q(idxT:idyT)));
            else
                % Extract all properties at the core PV depth level
                pv_data(i,1)=CT(q==min(q(idxT:idyT)));
                pv_data(i,2)=SA(q==min(q(idxT:idyT)));
                pv_data(i,3)=p_dens(q==min(q(idxT:idyT)));
                pv_data(i,4)=CT_grad(q==min(q(idxT:idyT)));
                pv_data(i,5)=p_dens_grad(q==min(q(idxT:idyT)));
                %                 core_aou(i)=tmp_aou(q==min(q(idxT:idyT)));
            end
            
        catch % If there is no STMW layer this is triggered
            low_pv(i,1:3)=NaN;
            low_bv(i,1:2)=NaN;
            display(['Low_PV ',num2str(i),' failed']);
        end
    else
        low_pv(i,1:3)=NaN;
        low_bv(i,1:2)=NaN;
        pv_data(i,:)=NaN;
        no_layer=no_layer+1;
    end
    
    % Create pv contour
    try
        pv_cont(:,i)=interp1(pres,q,5:5:600);
    catch
        pv_cont(:,i)=NaN;
    end
    
    % Calculate MW intensity Qiu 2006
    if sum(s_msk)>2
        worthington_intensity(i,1)=...
            trapz(pres(idx:idy),(inpaint_nans(2e-10-q(idx:idy))));
    end
    % Create temperature gradient per 100m
    if sum(s_msk)>2
        w_grad(i,1)=((max(CT(s_msk))-min(CT(s_msk)))...
            /stevens_thickness(i))*100;
    end

end


% idx=stevens_thickness==0;
% mtime(idx)=[];
% pv_data(idx,:)=[];
% lat(idx)=[];
% lon(idx)=[];
% low_pv(idx,:)=[];
% worthington_intensity(idx)=[];
% stevens_thickness(idx,:)=[];
month=unique(time);
clearvars -except month pv_data low_pv worthington_intensity stevens_thickness

save('EN4g10_BATS.mat');

%% ALL stmw data
clear
load EN4g10_data_V2.mat

depth=[0:2:600]';

lat_list=unique(lat_vec);
lon_list=unique(lon_vec);

[lon_grid,lat_grid]=meshgrid(lon_list,lat_list);
thick_grid=NaN([size(lon_grid) length(1955:.25:2019)]);
ulim_grid=NaN([size(lon_grid) length(1955:.25:2019)]);
llim_grid=NaN([size(lon_grid) length(1955:.25:2019)]);
PV_grid=NaN([size(lon_grid) length(1955:.25:2019)]);
intensity_grid=NaN([size(lon_grid) length(1955:.25:2019)]);
stemp_grid=NaN([size(lon_grid) length(1955:.25:2019)]);
hc_grid=NaN([size(lon_grid) length(1958:.25:2018)]);
wtemp_grid=NaN([size(lon_grid) length(1955:2019)]);
outcrop_grid=zeros([size(lon_grid) length(1955:2019)]);

tic
% Find how many profiles have EDW outcropping and run tests to ensure they
% meet EDW criteria
for i=1:length(lon_list)
    for ii=1:length(lat_list)
        
        idx=lon_vec==lon_list(i) & lat_vec==lat_list(ii);
        tmp_time=time(idx);
        tmp_temp=temp_grid(:,idx);
        tmp_salt=salt_grid(:,idx);
        pres=gsw_p_from_z(depth*-1,lat_list(ii));
        count=0;
        count2=0;

        for iii=1955.5:.25:2019
            i
            ii
            count
            
%             if lat_list(ii)>32 && lat_list(ii)<35 && ...
%                     lon_list(i)>-64 && lon_list(i)<-60 && iii==1967
%                 keyboard
%             end
%             
            count=count+1;
            time_idx=tmp_time>iii & tmp_time<iii+.25;
            mn_temp=mean(tmp_temp(:,time_idx),2,'omitnan')-273.15;
            mn_salt=mean(tmp_salt(:,time_idx),2,'omitnan');
            SA=gsw_SA_from_SP(mn_salt,pres,lat_list(ii),lon_list(i));
            CT=gsw_CT_from_pt(SA,mn_temp);
            
            % Store surface (<15m) temperature values
            stemp_grid(ii,i,count)=nanmean(CT(depth<15));
            
            % Calculate surface (10-150m ) heat content
            idxd=depth>=10 & depth<150;
            hc_grid(ii,i,count)=trapz(depth(idxd),gsw_enthalpy(SA(idxd),...
            CT(idxd),depth(idxd)));
        
            % Calculate density from variables
            p_dens=gsw_sigma0(SA,CT);
            p_dens=inpaint_nans(p_dens);
            
            % Calculate Worthington thickness (thickness between 17 and 19
            % degree isotherms)
            [~,idy]=min(abs(17-CT));
            [~,idx]=min(abs(19-CT));
            
            % Calculate "Stevens thickness" (thickness between 17 and 19
            % degree isotherms and PV < 1x10^-10)
            [~,q,p_ave]=sw_bfrq(mn_salt,gsw_t_from_CT(SA,CT,pres),pres,lat_list(ii));
            q=interp1(p_ave,q,pres);
            
            % Filter PV signal with a moving mean to remove the high frequency variability
            q=movmean(q,35);
            msk=zeros(length(CT),1);
            msk(idx:idy)=1;
            s_msk=logical(msk.*(q<1e-10));
            
            if sum(s_msk)>2
                thick_grid(ii,i,count)=max(depth(s_msk))-...
                    min(depth(s_msk));
                ulim_grid(ii,i,count)=min(depth(s_msk));
                llim_grid(ii,i,count)=max(depth(s_msk));
            else
                thick_grid(ii,i,count)=0;
            end
            
            % Calulate core pv between 16oC isotherms
            if sum(s_msk)>2
                try
                    [~,idyT]=min(abs(16-CT));[~,idxT]=min(abs(20-CT));
                    PV_grid(ii,i,count)=min(q(idxT:idyT));
                    
                catch % If there is no STMW layer this is triggered
                    disp(['Low_PV ',num2str(i),' failed']);
                end
            end
            
            % Calculate MW intensity Qiu 2006
            if sum(s_msk)>2
                intensity_grid(ii,i,count)=...
                    trapz(pres(idx:idy),(inpaint_nans(2e-10-q(idx:idy))));
            end
            
            if rem(iii,1)==0
                count2=count2+1;
                
                time_idx=tmp_time>iii+0.0847 & tmp_time<iii+0.3306;
                mn_temp=mean(tmp_temp(:,time_idx),2,'omitnan')-273.15;
                mn_salt=mean(tmp_salt(:,time_idx),2,'omitnan');
                SA=gsw_SA_from_SP(mn_salt,pres,lat_list(ii),lon_list(i));
                CT=gsw_CT_from_pt(SA,mn_temp);
                
                % Store Feb-Mar surface (<15m) temperature values
                wtemp_grid(ii,i,count2)=nanmean(CT(depth<15));
                              
%                 if lon_list(i)>-63 && lon_list(i)<-65 && lat_list(ii)>31 && latlist(ii)<33
%                     keyboard
%                 end
                
                if wtemp_grid(ii,i,count2)>=17 && wtemp_grid(ii,i,count2)<=19 && sum(s_msk)>2
                    outcrop_grid(ii,i,count2)=thick_grid(ii,i,count2)*...
                        gsw_distance([lon_list(i)-0.5 lon_list(i)+0.5],...
                        [lat_list(ii) lat_list(ii)])*gsw_distance([lon_list(i) lon_list(i)],...
                        [lat_list(ii)-0.5 lat_list(ii)+0.5]);
                end
            end
        end
    end
end

time=1955:.25:2019;
clearvars -except years PV_grid thick_grid ulim_grid llim_grid stemp_grid intensity_grid lat_grid lon_grid wtemp_grid outcrop_grid hc_grid
save('EN4g10_all_V3.mat');

toc

