%% Read ARGO data
addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))

cd C:\Users\samst\Documents\BATS_VPN\SWS_STMW\Code\argo\

% east all
dir_struct1=dir('C:\Users\samst\Documents\BATS_VPN\SWS_STMW\Code\argo\data\*argo*'); 
% west
% dir_struct2=dir('C:\Users\samst\Documents\BATS_VPN\ARGO\raw_profiles\ALL\*argo*');

temp_grid=zeros(length(0:2:600),20000);
salt_grid=zeros(length(0:2:600),20000);
lat_vec=zeros(20000,1);
lon_vec=zeros(20000,1);
time_vec=zeros(20000,1);
count=1;
error=0;

for i=1:length(dir_struct1)
    
    ds=dir_struct1(i);
       
    try
        temp=ncread(ds.name,'TEMP_ADJUSTED');
        salt=ncread(ds.name,'PSAL_ADJUSTED');
        pres=ncread(ds.name,'PRES_ADJUSTED');
        lat=ncread(ds.name,'LATITUDE');
        lon=ncread(ds.name,'LONGITUDE');
        time=ncread(ds.name,'JULD');
        
        for ii=1:size(temp,2)
            msk=~isnan(pres(:,ii));
            temp_grid(:,count)=interp1(pres(msk,ii),temp(msk,ii),0:2:600);
            salt_grid(:,count)=interp1(pres(msk,ii),salt(msk,ii),0:2:600);
            lat_vec(count)=lat(ii);
            lon_vec(count)=lon(ii);
            time_vec(count)=time(ii)+datenum(1950,01,01);
            count=count+1;
        end
        
    catch
        disp('profile failed');
        error=error+1;
    end
end

zmsk=time_vec==0;
time_vec(zmsk)=[];
lat_vec(zmsk)=[];
lon_vec(zmsk)=[];
temp_grid(:,zmsk)=[];
salt_grid(:,zmsk)=[];

[mtime,idx]=sort(time_vec);
lat=lat_vec(idx);
lon=lon_vec(idx);
temp=temp_grid(:,idx);
salt=salt_grid(:,idx);

clearvars -except mtime lat lon temp salt error 
save('C:\Users\samst\Documents\BATS_VPN\SWS_STMW\Code\argo\data\argo_profiles.mat');

%% STMW Analysis
clear
load C:\Users\samst\Documents\BATS_VPN\SWS_STMW\Code\argo\data\argo_profiles.mat

pres=[0:2:600]';
worthington_intensity=zeros(size(temp,2),1);
pv_data=zeros(size(temp,2),5);
w_grad=zeros(size(temp,2),1);
mld=zeros(size(temp,2),3);
temp_cont=zeros(600/5,size(temp,2));
sal_cont=zeros(600/5,size(temp,2));
dens_cont=zeros(600/5,size(temp,2));
pv_cont=zeros(600/5,size(temp,2));
tgrad_cont=zeros(600/5,size(temp,2));
heat_content=zeros(size(temp,2),1);
stmw_heat_content=zeros(size(temp,2),1);
stmw_temp=zeros(size(temp,2),1);
stevens_thickness=zeros(size(temp,2),1);
stevens_limits=zeros(size(temp,2),2);
outcropping=zeros(size(temp,2),1);
core_aou=zeros(size(temp,2),1);
time=zeros(size(temp,2),1);
no_layer=0;
error=[];

% Find how many profiles have EDW outcropping and run tests to ensure they
% meet EDW criteria

for i=1:size(temp,2)
    
    % Bad data ids
    try        
        display(i);
        
        depth=gsw_z_from_p(pres,lat(i))*-1;
        
        % Calculate TEOS-10 variables
        SA=gsw_SA_from_SP(salt(:,i),pres,lon(i),lat(i));
        CT=gsw_CT_from_t(SA,temp(:,i),pres);
        
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
        [~,q,p_ave]=sw_bfrq(salt(:,i),temp(:,i),pres,lat(i));
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
    catch
        error=[error;i];
    end
    
    time(i)=cnv_date(mtime(i));
end

idx=stevens_thickness==0;
mtime(idx)=[];
pv_data(idx,:)=[];
lat(idx)=[];
lon(idx)=[];
low_pv(idx,:)=[];
worthington_intensity(idx)=[];
stevens_thickness(idx,:)=[];
