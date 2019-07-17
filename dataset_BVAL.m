% Create dataset for study
addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))

internal=load('CTD_internal.mat');
CTD_data=internal.CTD_data;
clearvars -except CTD_data

CTD_data(floor(CTD_data(:,1)/10000000)~=5,:)=[];
CTD_data=sortrows(CTD_data,2);

% Create castlist
[d,I]=unique(CTD_data(:,1),'first');
castlist=CTD_data(sort(I),1);

% Preallocate arrays
time=zeros(length(castlist),1);
worthington=zeros(length(castlist),1);
worthington_limits=zeros(length(castlist),2);
worthington_intensity=zeros(length(castlist),1);
pv_data=zeros(length(castlist),5);
ed=zeros(length(castlist),1);
ed_limits=zeros(length(castlist),2);
w_grad=zeros(length(castlist),1);
mld=zeros(length(castlist),3);
temp_cont=zeros(600/5,length(castlist));
sal_cont=zeros(600/5,length(castlist));
dens_cont=zeros(600/5,length(castlist));
pv_cont=zeros(600/5,length(castlist));
heat_content=zeros(length(castlist),1);
stmw_heat_content=zeros(length(castlist),1);
stmw_temp=zeros(length(castlist),1);
stevens_thickness=zeros(length(castlist),1);
stevens_limits=zeros(length(castlist),2);
outcropping=zeros(length(castlist),1);
core_aou=zeros(length(castlist),1);
coords=zeros(length(castlist),2);

tic
for i=1:length(castlist)
    
    %find individual cast
    rel_cast=CTD_data((CTD_data(:,1)==castlist(i)),:);
    
    if max(rel_cast(:,5))>490 && min(rel_cast(:,5))<100
        display(rel_cast(1,1));
        
        % Some casts are not depth sorted for some reason or have repeated
        % depth bins... make sure they are/don't
        [~,idx]=unique(rel_cast(:,5));
        rel_cast=rel_cast(idx,:);
        
         % Log coordinates
        coords(i,:)=[rel_cast(1,4)*-1 rel_cast(1,3)];
        
        %create time vector
        time(i)=rel_cast(1,2);
        
        % Calculate Coriolis parameter
        f=gsw_f(rel_cast(1,3));
        
        % Calculate necessary variables
        SA=gsw_SA_from_SP(rel_cast(:,9),rel_cast(:,5),rel_cast(1,4),rel_cast(1,3));
        CT=gsw_CT_from_t(SA,rel_cast(:,7),rel_cast(:,5));
        
        % Calculate density from variables
        p_dens=gsw_sigma0(SA,CT);
        p_dens=inpaint_nans(p_dens);
        
        % Calculate Worthington thickness (thickness between 17 and 19
        % degree isotherms)
        [~,idy]=min(abs(17-CT));
        mw_low=rel_cast(idy,6);
        
        worthington_limits(i,2)=mw_low;
        
        [~,idx]=min(abs(19-CT));
        mw_high=rel_cast(idx,6);
        worthington_limits(i,1)=mw_high;
        worthington(i,1)=mw_low-mw_high;
        
        % Calculate "Stevens thickness" (thickness between 17 and 19
        % degree isotherms and PV < 1x10^-10). Filter PV signal with a
        % moving mean to remove the high frequency variability
        [~,q,p_ave]=sw_bfrq(rel_cast(:,9),rel_cast(:,7),rel_cast(:,5),rel_cast(1,3));        
        q=interp1(p_ave,q,rel_cast(:,5));
        q=movmean(q,35);
        msk=zeros(length(CT),1);
        msk(idx:idy)=1;
        s_msk=logical(msk.*(q<1e-10));
        
        msk=zeros(length(CT),1);
        msk(idx:idy)=1;
        s_msk=logical(msk.*(q<1e-10));
        
        if sum(s_msk)>2
            stevens_thickness(i)=max(rel_cast((s_msk),6))-...
                min(rel_cast((s_msk),6));
            stevens_limits(i,:)=[min(rel_cast(s_msk,6)) ...
                max(rel_cast((s_msk),6))];
        end
        
%         % QC check: uncomment to plot
%         figure
%         subplot(1,2,1);
%         plot(CT,rel_cast(:,5));
%         set(gca,'ydir','reverse');
%         hold on
%         plot(CT(s_msk),rel_cast(s_msk,5),'linewidth',2);
%         plot(CT(idx:idy),rel_cast(idx:idy,5),'--','linewidth',2);
%         ylim([0 600]);
%         xlim([15 22]);
%         xlabel('Temp./ \circ C');
%         subplot(1,2,2);
%         plot(q,rel_cast(:,5));
%         set(gca,'ydir','reverse');
%         hold on
%         plot(q(s_msk),rel_cast(s_msk,5),'linewidth',2);
%         plot(q(idx:idy),rel_cast(idx:idy,5),'--','linewidth',2);
%         ylim([0 600]);
%         min(q(idx:idy))
%         legend('smth prof.','19-17C isotherms','Stevens thickness',...
%             'location','best');
%         xlabel('PV/ m^{-1} s^{-1}');
%         pause
%         close all
        
        % Calulate MLD
        try
            mld(i,1)=ra_mld(SA,CT,rel_cast(:,6),0.5); % 1° MLD temp as per Hanawa and Hoshino (1988) and Sufa and Hanawa (1990)
            mld(i,2)=ra_mld(SA,CT,rel_cast(:,6),1);
            mld(i,3)=rel_cast(1,1);
        catch
            mld(i,:)=NaN;
        end
        
        % Create temperature contour
        try
            temp_cont(:,i)=interp1(rel_cast(:,6),CT,5:5:600);
        catch
            temp_cont(:,i)=NaN;
        end
        
        % Create salinity contour
        try
            sal_cont(:,i)=interp1(rel_cast(:,6),SA,5:5:600);
        catch
            sal_cont(:,i)=NaN;
        end
        
        % Calculate enthalpy
        idxh=~isnan(CT) & ~isnan(SA) &...
            ~isnan(rel_cast(:,5)) & rel_cast(:,5)<150 &...
            rel_cast(:,5)>10;
        
        
        % Create surface integrated heat content <150m
        heat_content(i)=trapz(rel_cast(idxh,5),gsw_enthalpy_CT_exact(SA(idxh),...
            CT(idxh),rel_cast(idxh,5)));
        
        % Create STMW integrated heat content
        if sum(s_msk)>2
            stmw_heat_content(i)=trapz(rel_cast(s_msk,5),gsw_enthalpy_CT_exact(SA(s_msk),...
                CT(s_msk),rel_cast(s_msk,5)));
            
            stmw_temp(i)=nanmean(CT(s_msk));
        else
            stmw_heat_content(i)=NaN;
            
            stmw_temp(i)=NaN;
        end
        
        % Outcropping?
        OC_idx=rel_cast(:,5)<10;
        if sum(OC_idx)>2 && nanmean(CT(OC_idx))<19
            outcropping(i)=1;
        end
        
        % Create temperature and density gradients
        CT_grad=diff(CT);
        p_dens_grad=diff(p_dens);
        
        % mean AOU in the STMW zone
        PT=gsw_pt_from_CT(SA,CT);
        tmp_aou=aou_new(SA,PT,rel_cast(:,10));
        
        %         [~,idy_c]=min(abs(17.5-CT));
        %         [~,idx_c]=min(abs(18.5-CT));
        idx_c=idx;
        idy_c=idy;
        
        % Calulate core pv
        try
            low_pv(i,1)=min(q(idx_c:idy_c)); low_pv(i,2)=rel_cast(q==min(q(idx_c:idy_c)),6); low_pv(i,3)=rel_cast(1,1);
            
            % Interpolate all values to the core PV depth level
            pv_data(i,1)=CT(q==min(q(idx_c:idy_c)));
            pv_data(i,2)=SA(q==min(q(idx_c:idy_c)));
            pv_data(i,3)=p_dens(q==min(q(idx_c:idy_c)));
            pv_data(i,4)=CT_grad(q==min(q(idx_c:idy_c)));
            pv_data(i,5)=p_dens_grad(q==min(q(idx_c:idy_c)));
            core_aou(i)=tmp_aou(q==min(q(idx_c:idy_c)));
        catch
            low_pv(i,1:3)=NaN;
            low_bv(i,1:2)=NaN;
            display(['Low_PV ',num2str(rel_cast(1,1)),' failed']);
        end
        
        % Create pv contour
        try
            pv_cont(:,i)=interp1(rel_cast(:,6),q,5:5:600);
        catch
            pv_cont(:,i)=NaN;
        end
        
         % Calculate MW intensity Qiu 2006
        worthington_intensity(i,1)=...
            trapz(rel_cast(idx:idy,5),(inpaint_nans(2e-10-q(idx:idy))));
        
        
        % Create temperature gradient per 100m
        if sum(s_msk)>2
            w_grad(i,1)=((max(CT(s_msk))-min(CT(s_msk)))...
                /stevens_thickness(i))*100;
        end
        
        % Create ebbesmeyer density thickness
        ed_limits(i,1)=interp1(p_dens(~isnan(p_dens)),...
            rel_cast(~isnan(p_dens),6),26.33);
        ed_limits(i,2)=interp1(p_dens(~isnan(p_dens)),...
            rel_cast(~isnan(p_dens),6),26.51);
        ed(i,1)=ed_limits(i,2)-ed_limits(i,1);
    end
    
end
toc

% Remove all datapoints that did not pass original conditions
low_pv_final=low_pv((time>0),:);
pv_data_final=pv_data((time>0),:);
time_final=time(time>0);
w_grad_final=w_grad((time>0),:);
mld_final=mld((time>0),:);
temp_cont_final=temp_cont(:,(time>0));
pv_cont_final=pv_cont(:,(time>0));
worthington_intensity_final=worthington_intensity((time>0),:);
heat_content_final=heat_content((time>0),:);
stmw_heat_content_final=stmw_heat_content((time>0),:);
stevens_thickness_final=stevens_thickness((time>0),:);
outcropping_final=outcropping((time>0),:);
stevens_limits_final=stevens_limits((time>0),:);
worthington_final=worthington((time>0),:);
aou_final=core_aou((time>0),:);
ed_final=ed((time>0),:);
ed_limits_final=ed_limits((time>0),:);
coords_final=coords((time>0),:);
