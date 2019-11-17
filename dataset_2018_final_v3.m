% Create dataset for study
addpath(genpath('C:\Users\samst\UBC'))
addpath(genpath('C:\Users\samst\Documents'))
% internal=load('CTD_internal.mat');
internal=load('C:\Users\samst\Documents\BATS_VPN\CTD_internal_oct2019.mat');
CTD_data=internal.CTD_data;
clearvars -except CTD_data

% Remove non HS and BATS data and problematic cruises, sort data by date
CTD_data((CTD_data(:,1)>30000000 & CTD_data(:,1)<60000000),:)=[];
CTD_data((CTD_data(:,1)>70000000),:)=[];
CTD_data((CTD_data(:,1)==61324002),:)=[];
CTD_data((CTD_data(:,1)==10215020),:)=[];
CTD_data((CTD_data(:,1)==61203001),:)=[];
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
gmax=zeros(length(castlist),1);
mld=zeros(length(castlist),3);
temp_cont=zeros(600/5,length(castlist));
sal_cont=zeros(600/5,length(castlist));
dens_cont=zeros(600/5,length(castlist));
dens_grad_cont=zeros(600/5,length(castlist));
pv_cont=zeros(600/5,length(castlist));
tgrad_cont=zeros(600/5,length(castlist));
heat_content=zeros(length(castlist),1);
stmw_heat_content=zeros(length(castlist),1);
stmw_temp=zeros(length(castlist),1);
stevens_thickness=zeros(length(castlist),1);
stevens_limits=zeros(length(castlist),2);
outcropping=zeros(length(castlist),1);
core_aou=zeros(length(castlist),1);
no_layer=0;
s_spice=zeros(length(castlist),1);
UO_heat_content=zeros(length(castlist),1);
thickness_def=zeros(length(castlist),5);
pv_def=zeros(length(castlist),5);
t_def=NaN(length(castlist),5);

tic
for i=1:length(castlist)
    
    %find individual cast
    rel_cast=CTD_data((CTD_data(:,1)==castlist(i)),:);
    
    if max(rel_cast(:,5))>490 && min(rel_cast(:,5))<100 &&...
        sum(isnan(rel_cast(:,7)))<5 && sum(isnan(rel_cast(:,9)))<5 %&&...
        %rel_cast(1,2)<=2018.489785007610+.000000000001
    
        display(rel_cast(1,1));
                
        % Some casts are not depth sorted for some reason or have repeated
        % depth bins... make sure they are/don't
        [~,idx]=unique(rel_cast(:,5));
        rel_cast=rel_cast(idx,:);
        
        %create time vector
        time(i)=rel_cast(1,2);
        
        % Calculate Coriolis parameter
        f=gsw_f(rel_cast(1,3));
        
        % Calculate TEOS-10 variables
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
        % degree isotherms and PV < 1x10^-10)
        [~,q,p_ave]=sw_bfrq(rel_cast(:,9),rel_cast(:,7),rel_cast(:,5),rel_cast(1,3));
        q=interp1(p_ave,q,rel_cast(:,5));

        % Filter PV signal with a moving mean to remove the high frequency variability
        q=movmean(q,35);
        msk=zeros(length(CT),1);
        msk(idx:idy)=1;
        s_msk=logical(msk.*(q<1e-10));
        
        if sum(s_msk)>2
            stevens_thickness(i)=max(rel_cast((s_msk),6))-...
                min(rel_cast((s_msk),6));
            stevens_limits(i,:)=[min(rel_cast(s_msk,6)) ...
                max(rel_cast((s_msk),6))];
        end
        
        % Calculate spice
        tmp_spice=spice(rel_cast(:,5),CT,SA);
                
        % QC check plot (uncomment if needed)
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
        
        % Create temperature and density gradients
        CT_grad=diff(CT);
        p_dens_grad=diff(p_dens);
        
          % Create temp. gradient contour
        try
            tgrad_cont(:,i)=interp1(p_ave,CT_grad,5:5:600);
        catch
            tgrad_cont(:,i)=NaN;
        end

        % Create dens. gradient contour
        try
            d_grad_cont(:,i)=interp1(p_ave,p_dens_grad,5:5:600);
        catch
            d_grad_cont(:,i)=NaN;
        end
        
        % Calculate enthalpy
        idxh=~isnan(CT) & ~isnan(SA) &...
            ~isnan(rel_cast(:,5)) & rel_cast(:,5)<150 &...
            rel_cast(:,5)>10;
        
        % Create surface integrated heat content <150m
        heat_content(i)=trapz(rel_cast(idxh,5),gsw_enthalpy(SA(idxh),...
            CT(idxh),rel_cast(idxh,5)));
        
        idxd=~isnan(CT) & ~isnan(SA) &...
            ~isnan(rel_cast(:,5)) & rel_cast(:,5)<700 &...
            rel_cast(:,5)>10;
        
        % Create upper ocean integrated heat content <700m
        UO_heat_content(i)=trapz(rel_cast(idxd,5),gsw_enthalpy(SA(idxd),...
            CT(idxd),rel_cast(idxd,5)));
        
        % Create STMW integrated heat content
        if sum(s_msk)>2
            stmw_heat_content(i)=trapz(rel_cast(s_msk,5),gsw_enthalpy(SA(s_msk),...
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
        
%         % mean AOU in the STMW zone
%         PT=gsw_pt_from_CT(SA,CT);
%         tmp_aou=aou_new(SA,PT,rel_cast(:,10));
        
        % Calulate core pv between 16oC isotherms
        if sum(s_msk)>2
            try
                [~,idyT]=min(abs(16-CT));[~,idxT]=min(abs(20-CT));
                low_pv(i,1)=min(q(idxT:idyT)); low_pv(i,2)=rel_cast(q==min(q(idxT:idyT)),6); low_pv(i,3)=rel_cast(1,1);
                
                if low_pv(i,1)<0
                    low_pv(i,1)=min(q(idx:idy)); low_pv(i,2)=rel_cast(q==min(q(idx:idy)),6); low_pv(i,3)=rel_cast(1,1);
                    % Extract all properties at the core PV depth level
                    pv_data(i,1)=CT(q==min(q(idx:idy)));
                    pv_data(i,2)=SA(q==min(q(idx:idy)));
                    pv_data(i,3)=p_dens(q==min(q(idx:idy)));
                    pv_data(i,4)=C_grad(q==min(q(idx:idy)));
                    pv_data(i,5)=p_dens_grad(q==min(q(idx:idy)));
                    s_spice(i)=tmp_spice(q==min(q(idx:idy)));
                    keyboard
                    %                 core_aou(i)=tmp_aou(q==min(q(idxT:idyT)));
                else
                    % Extract all properties at the core PV depth level
                    pv_data(i,1)=CT(q==min(q(idxT:idyT)));
                    pv_data(i,2)=SA(q==min(q(idxT:idyT)));
                    pv_data(i,3)=p_dens(q==min(q(idxT:idyT)));
                    pv_data(i,4)=CT_grad(q==min(q(idxT:idyT)));
                    pv_data(i,5)=p_dens_grad(q==min(q(idxT:idyT)));
                    s_spice(i)=tmp_spice(q==min(q(idxT:idyT)));
                    %                 core_aou(i)=tmp_aou(q==min(q(idxT:idyT)));
                end
                
            catch % If there is no STMW layer this is triggered
                low_pv(i,1:3)=NaN;
                low_bv(i,1:2)=NaN;
                s_spice(i)=NaN;
                display(['Low_PV ',num2str(rel_cast(1,1)),' failed']);
            end
        else
            low_pv(i,1:3)=NaN;
            low_bv(i,1:2)=NaN;
            pv_data(i,:)=NaN;
            s_spice(i)=NaN;
            no_layer=no_layer+1;
        end
        
        % Create pv contour
        try
            pv_cont(:,i)=interp1(rel_cast(:,5),q,5:5:600);
        catch
            pv_cont(:,i)=NaN;
        end
        
        % Create dens contour
        try
            dens_cont(:,i)=interp1(rel_cast(:,5),p_dens,5:5:600);
        catch
            dens_cont(:,i)=NaN;
        end

        % Calculate MW intensity Qiu 2006
        worthington_intensity(i,1)=...
            trapz(rel_cast(idx:idy,5),(inpaint_nans(2e-10-q(idx:idy))));
        
        % Create temperature gradient per 100m
        if sum(s_msk)>2
            w_grad(i,1)=((max(CT(s_msk))-min(CT(s_msk)))...
                /stevens_thickness(i))*100;
        end

        % find Gmax (maximum density stratification in surface
        try
            [gmax(i),tmp]=max(smooth(p_dens_grad(rel_cast(1:end-1,5)<200),5));
        catch
            gmax(i)=NaN;
        end
        
        % Create ebbesmeyer density thickness
        ed_limits(i,1)=interp1(p_dens(~isnan(p_dens)),...
            rel_cast(~isnan(p_dens),6),26.33);
        ed_limits(i,2)=interp1(p_dens(~isnan(p_dens)),...
            rel_cast(~isnan(p_dens),6),26.51);
        ed(i,1)=ed_limits(i,2)-ed_limits(i,1);
        
        % Calculate thickness based on less stringent range (thickness between 17 and 19
        % degree isotherms, PV <1.5 x 10^-10)
        [~,idy]=min(abs(16-CT));
        [~,idx]=min(abs(20-CT));
        msk=zeros(length(CT),1);
        msk(idx:idy)=1;
        s_msk=logical(msk.*(q<2e-10));
        
        if sum(s_msk)>2
            thickness_def(i,1)=max(rel_cast((s_msk),6))-...
                min(rel_cast((s_msk),6));
             try
                pv_def(i,1)=min(q(s_msk));
                t_def(i,1)=CT(q==min(q(s_msk)));
             end
            
        end
        
        % Calculate thickness with density limit (thickness between 17 and 19
        % degree isotherms, PV <1.5 x 10^-10, density 26.2-26.6 isopycnals)
        [~,idy]=min(abs(17.5-CT));
        [~,idx]=min(abs(18.5-CT));
        msk=zeros(length(CT),1);
        msk(idx:idy)=1;
        s_msk=logical(msk.*(q<1e-10));
        
        if sum(s_msk)>2
            thickness_def(i,2)=max(rel_cast((s_msk),6))-...
                min(rel_cast((s_msk),6));
            try
                pv_def(i,2)=min(q(s_msk));
                t_def(i,2)=CT(q==min(q(s_msk)));
             end
        end
        
        % Calculate thickness density 26.3-26.5 isopycnals)
        s_msk=p_dens>26.3 & p_dens<26.5;
        
        if sum(s_msk)>2
            thickness_def(i,3)=max(rel_cast((s_msk),6))-...
                min(rel_cast((s_msk),6));
            try
                pv_def(i,3)=min(q(s_msk));
                t_def(i,3)=CT(q==min(q(s_msk)));
            end
        end
        
        % Calculate thickness with isotherm, density and PV limit (PV <1 x
        % 10^-10, density 26.3-26.5 isopycnals, 17-19 isotherm)
        [~,idy]=min(abs(17-CT));
        [~,idx]=min(abs(19-CT));
        msk=zeros(length(CT),1);
        msk(idx:idy)=1;
        s_msk=logical(msk.*(q<1e-10).*(p_dens>26.3 & p_dens<26.5));
                
        if sum(s_msk)>2
            thickness_def(i,4)=max(rel_cast((s_msk),6))-...
                min(rel_cast((s_msk),6));
            try
                pv_def(i,4)=min(q(s_msk));
                t_def(i,4)=CT(q==min(q(s_msk)));
            end
        end
        
        % Calculate thickness with isotherm, density and PV limit (PV <1 x
        % 10^-10, density 26.3-26.5 isopycnals, 17-19 isotherm)
        [~,idy]=min(abs(17-CT));
        [~,idx]=min(abs(19-CT));
        msk=zeros(length(CT),1);
        msk(idx:idy)=1;
        CTgi=interp1(p_ave,CT_grad,rel_cast(:,5))/2;
        s_msk=logical(msk.*(smooth(CTgi/2,10)<0.06).*(p_dens>26.2 & p_dens<26.7));
                
        if sum(s_msk)>2
            thickness_def(i,5)=max(rel_cast((s_msk),6))-...
                min(rel_cast((s_msk),6));
            try
                pv_def(i,5)=min(q(s_msk));
                t_def(i,5)=CT(q==min(q(s_msk)));
            end
        end
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
dens_cont_final=dens_cont(:,(time>0));
worthington_intensity_final=worthington_intensity((time>0),:);
heat_content_final=heat_content((time>0),:);
stmw_heat_content_final=stmw_heat_content((time>0),:);
stevens_thickness_final=stevens_thickness((time>0),:);
outcropping_final=outcropping((time>0),:);
stevens_limits_final=stevens_limits((time>0),:);
worthington_final=worthington((time>0),:);
% aou_final=core_aou((time>0),:);
ed_final=ed((time>0),:);
ed_limits_final=ed_limits((time>0),:);
tgrad_cont_final=tgrad_cont(:,(time>0));
d_grad_cont_final=d_grad_cont(:,(time>0));
s_spice_final=s_spice(time>0);
UO_heat_content_final=UO_heat_content(time>0);
gmax_final=gmax(time>0);
thickness_def_final=thickness_def(time>0,:);
pv_def_final=pv_def(time>0,:);
t_def_final=t_def(time>0,:);
