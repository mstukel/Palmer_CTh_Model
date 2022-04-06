function [norm_misfit] = PhyDetFPmodelFirstOrder(k_sorp_phy,k_sorp_det,k_desorp,Eg_Th,starttimes,endtimes,toplot)

%Eg_Th = Fraction of thorium consumed by krill that passes to their fecal pellets

% clearvars
% close all
% k_sorp_phy = 1/30;
% k_sorp_det = 1/30;
% k_desorp = 1/365;
% Eg_Th = 0.5;
% starttimes = [9.5,50];
% endtimes = [40,137];
% toplot = 1;

Eg_C = 0.3;  %Fraction of carbon consumed by prill that passes to their fecal pellets (i.e., egestion efficiency = 1 - assimilation efficiency)

%Diagnostics
CompileValidationData=0;

%Setting up boundaries and forcing
load('Phy-Det-FP model inputs 1.mat')
load('Phy-Det-FP model inputs 2.mat')
load('Phy-Det-FP model inputs 3.mat')
Phy_hold = Phy(1:find(z==65),:);
Det_hold = Det(1:find(z==65),:);
PTh_hold = PTh(1:find(z==65),:);
PhyTh_hold = PTh_hold .* Phy_hold./(Phy_hold+Det_hold);
DetTh_hold = PTh_hold .* Det_hold./(Phy_hold+Det_hold);
dTh_hold = Th(1:find(z==65),:)-PTh_hold;
dt = time(2)-time(1);
C14PP = C14PP(1:find(z==65),:);
Grazing = Grazing(1:find(z==65),:);
MortPhy2Part = MortPhy2Part(1:find(z==65),:);
MortPhy2Diss = MortPhy2Diss(1:find(z==65),:);
Remin = Remin(1:find(z==65),:);
omega = omega(1:find(z==65),:);
Kz = Kz(1:find(z_edge==67.5),:);
z=z(1:find(z==65));
z_edge=z_edge(1:find(z_edge==67.5));
z_thick=z_edge(2:end)-z_edge(1:end-1);

%Unchanging parameters
lam = 0.0287612937991679;


%Initializing variables
Phy_track = NaN(size(Phy_hold));
Det_track = NaN(size(Det_hold));
PhyTh_track = NaN(size(Phy_hold));
DetTh_track = NaN(size(Phy_hold));
dTh_track = NaN(size(Phy_hold));
POCSink_track = NaN(size(time));
ThSink_track = NaN(size(time));
FPCFlux_track = NaN(size(time));
FPThFlux_track = NaN(size(time));

%Model Run - First portion
Phy = Phy_hold(:,find(time==starttimes(1)));
Det = Det_hold(:,find(time==starttimes(1)));
PhyTh = PhyTh_hold(:,find(time==starttimes(1)));
DetTh = DetTh_hold(:,find(time==starttimes(1)));
dTh = dTh_hold(:,find(time==starttimes(1)));
clear Th_temp
for i=find(time==starttimes(1)):find(time==endtimes(1))

    %Sinking
    [Det_out,POCFlux] = sinking(Det,z_edge',omega(:,i),dt);
    Det = Det_out;
    
    [DetTh_out,PThFlux] = sinking(DetTh,z_edge',omega(:,i),dt);
    DetTh = DetTh_out;
    
    %Mixing
    [Coeff0,Coeff1,Coeff2,BottomCoeff]=CalculateMixingCoefficients(Det,z,z_edge,Kz(:,i),dt,Det_hold(end,i));
    [Det_out]=mixing_ftcs(Det,Det_hold(end,i),Coeff0,Coeff1,Coeff2);
    Det = Det_out;
    [Phy_out]=mixing_ftcs(Phy,Phy_hold(end,i),Coeff0,Coeff1,Coeff2);
    Phy = Phy_out;
    
    [PhyTh_out]=mixing_ftcs(PhyTh,PhyTh(end),Coeff0,Coeff1,Coeff2);
    PhyTh = PhyTh_out;
    [DetTh_out]=mixing_ftcs(DetTh,DetTh(end),Coeff0,Coeff1,Coeff2);
    DetTh = DetTh_out;
    
    [dTh_out]=mixing_ftcs(dTh,dTh_hold(end,i),Coeff0,Coeff1,Coeff2);
    dTh = dTh_out;

    %Primary Production
    Phy = Phy + C14PP(:,i)*dt;
    
    %Grazing
    FPCFlux_track(i) = sum(Grazing(:,i)*dt.*z_thick')*Eg_C;
    FPThFlux_track(i) = sum(Grazing(:,i)*dt.* PhyTh./Phy.*z_thick')*Eg_Th;
    dTh = dTh + Grazing(:,i)*dt*(1-Eg_Th) .* PhyTh./Phy;
    PhyTh = PhyTh - Grazing(:,i)*dt .* PhyTh./Phy;
    Phy = Phy - Grazing(:,i)*dt;
    
    %Phytoplankton Mortality
    dTh = dTh + MortPhy2Diss(:,i)*dt .* PhyTh./Phy;
    DetTh = DetTh + MortPhy2Part(:,i)*dt .* PhyTh./Phy;
    PhyTh = PhyTh - (MortPhy2Diss(:,i)+MortPhy2Part(:,i))*dt .* PhyTh./Phy;
    Phy = Phy - (MortPhy2Diss(:,i)+MortPhy2Part(:,i))*dt;
    Det = Det + MortPhy2Part(:,i)*dt;
    
    %Remineralization
    DetTh = DetTh - Remin(:,i)*dt .* DetTh./Det;
    dTh = dTh + Remin(:,i)*dt .* DetTh./Det;
    Det = Det - Remin(:,i)*dt;
    
    %Sorption
    PhyTh = PhyTh + k_sorp_phy*dTh*dt;
    dTh = dTh - k_sorp_phy*dTh*dt;
    DetTh = DetTh + k_sorp_det*dTh*dt;
    dTh = dTh - k_sorp_det*dTh*dt;
    
    %Desorption
    diff = k_desorp*PhyTh*dt;
    PhyTh = PhyTh - diff;
    dTh = dTh + diff;
    diff = k_desorp*DetTh*dt;
    DetTh = DetTh - diff;
    dTh = dTh + diff;
    
    %Th-234 Decay and production
    PhyTh = PhyTh - PhyTh*lam*dt;
    DetTh = DetTh - DetTh*lam*dt;
    dTh = dTh - dTh*lam*dt;
    dTh = dTh + U238(:,i)*lam*dt;
    
    Phy_track(:,i+1)=Phy;
    Det_track(:,i+1)=Det;
    PhyTh_track(:,i+1)=PhyTh;
    DetTh_track(:,i+1)=DetTh;
    dTh_track(:,i+1)=dTh;
    ThSink_track(i)=PThFlux(find(z==50));
    POCSink_track(i)=POCFlux(find(z==50));
end

%Model Run - Second portion
Phy = Phy_hold(:,find(time==starttimes(2)));
Det = Det_hold(:,find(time==starttimes(2)));
PhyTh = PhyTh_hold(:,find(time==starttimes(2)));
DetTh = DetTh_hold(:,find(time==starttimes(2)));
dTh = dTh_hold(:,find(time==starttimes(1)));
clear Th_temp
for i=find(time==starttimes(2)):find(time==endtimes(2))

    %Sinking
    [Det_out,POCFlux] = sinking(Det,z_edge',omega(:,i),dt);
    Det = Det_out;
    
    [DetTh_out,PThFlux] = sinking(DetTh,z_edge',omega(:,i),dt);
    DetTh = DetTh_out;
    
    %Mixing
    [Coeff0,Coeff1,Coeff2,BottomCoeff]=CalculateMixingCoefficients(Det,z,z_edge,Kz(:,i),dt,Det_hold(end,i));
    [Det_out]=mixing_ftcs(Det,Det_hold(end,i),Coeff0,Coeff1,Coeff2);
    Det = Det_out;
    [Phy_out]=mixing_ftcs(Phy,Phy_hold(end,i),Coeff0,Coeff1,Coeff2);
    Phy = Phy_out;
    
    [PhyTh_out]=mixing_ftcs(PhyTh,PhyTh(end),Coeff0,Coeff1,Coeff2);
    PhyTh = PhyTh_out;
    [DetTh_out]=mixing_ftcs(DetTh,DetTh(end),Coeff0,Coeff1,Coeff2);
    DetTh = DetTh_out;
    
    [dTh_out]=mixing_ftcs(dTh,dTh_hold(end,i),Coeff0,Coeff1,Coeff2);
    dTh = dTh_out;

    %Primary Production
    Phy = Phy + C14PP(:,i)*dt;
    
    %Grazing
    FPCFlux_track(i) = sum(Grazing(:,i)*dt.*z_thick')*Eg_C;
    FPThFlux_track(i) = sum(Grazing(:,i)*dt.* PhyTh./Phy.*z_thick')*Eg_Th;
    dTh = dTh + Grazing(:,i)*dt*(1-Eg_Th) .* PhyTh./Phy;
    PhyTh = PhyTh - Grazing(:,i)*dt .* PhyTh./Phy;
    Phy = Phy - Grazing(:,i)*dt;
    
    %Phytoplankton Mortality
    dTh = dTh + MortPhy2Diss(:,i)*dt .* PhyTh./Phy;
    DetTh = DetTh + MortPhy2Part(:,i)*dt .* PhyTh./Phy;
    PhyTh = PhyTh - (MortPhy2Diss(:,i)+MortPhy2Part(:,i))*dt .* PhyTh./Phy;
    Phy = Phy - (MortPhy2Diss(:,i)+MortPhy2Part(:,i))*dt;
    Det = Det + MortPhy2Part(:,i)*dt;
    
    %Remineralization
    DetTh = DetTh - Remin(:,i)*dt .* DetTh./Det;
    dTh = dTh + Remin(:,i)*dt .* DetTh./Det;
    Det = Det - Remin(:,i)*dt;
    
    %Sorption
    PhyTh = PhyTh + k_sorp_phy*dTh*dt;
    dTh = dTh - k_sorp_phy*dTh*dt;
    DetTh = DetTh + k_sorp_det*dTh*dt;
    dTh = dTh - k_sorp_det*dTh*dt;
    
    %Desorption
    diff = k_desorp*PhyTh*dt;
    PhyTh = PhyTh - diff;
    dTh = dTh + diff;
    diff = k_desorp*DetTh*dt;
    DetTh = DetTh - diff;
    dTh = dTh + diff;
    
    %Th-234 Decay and production
    PhyTh = PhyTh - PhyTh*lam*dt;
    DetTh = DetTh - DetTh*lam*dt;
    dTh = dTh - dTh*lam*dt;
    dTh = dTh + U238(:,i)*lam*dt;
    
    Phy_track(:,i+1)=Phy;
    Det_track(:,i+1)=Det;
    PhyTh_track(:,i+1)=PhyTh;
    DetTh_track(:,i+1)=DetTh;
    dTh_track(:,i+1)=dTh;
    ThSink_track(i)=PThFlux(find(z==50));
    POCSink_track(i)=POCFlux(find(z==50));
end

%Combining sinking POC and fecal pellets
SumSink_track = POCSink_track + FPCFlux_track;
SumThSink_track = ThSink_track + FPThFlux_track;

%Combining Phy and Det
POC_track = Phy_track + Det_track;
PTh_track = PhyTh_track + DetTh_track;

%Loading the field data for validation
if CompileValidationData==1
    load('..\..\..\MATLAB Synthesis Files\FieldData_ParticulateTh.mat')
    Field_PTh = table2array(FieldData_PTh);
    CTh = Field_PTh(:,5)/12./Field_PTh(:,3);  %C:Th umol/dpm
    CTh_err = CTh .* sqrt((Field_PTh(:,4)/100).^2 + (Field_PTh(:,7)./Field_PTh(:,5)).^2);
    for i=1:length(Field_PTh(:,1))
        ind = find(abs(time-Field_PTh(i,1))==min(abs(time-Field_PTh(i,1))));
        Field_PTh_Time_Indices(i,1) = ind(1);
        ind = find(abs(z-Field_PTh(i,2))==min(abs(z-Field_PTh(i,2))));
        Field_PTh_Depth_Indices(i,1) = ind(1);
    end
    Field_CTh=table(Field_PTh_Time_Indices,Field_PTh_Depth_Indices,CTh,CTh_err);
    Field_CTh=[FieldData_PTh(:,1:2),Field_CTh];
    load('..\..\..\MATLAB Synthesis Files\FieldData_Sedtrap.mat')
    Field_SedTrap = table2array(FieldData_SedTrap);
    for i=1:length(Field_SedTrap(:,1))
        ind = find(abs(time-Field_SedTrap(i,1))==min(abs(time-Field_SedTrap(i,1))));
        TimeStart_Indices(i,1) = ind(1);
        ind = find(abs(time-Field_SedTrap(i,2))==min(abs(time-Field_SedTrap(i,2))));
        TimeEnd_Indices(i,1) = ind(1);
    end
    temp=table(TimeStart_Indices,TimeEnd_Indices);
    Field_SedTrap=[temp,FieldData_SedTrap];
    origin = 'POCmodelFirstOrder.m';
    save('ValidationData','origin','Field_SedTrap','Field_CTh');
else
    load('ValidationData.mat')
end
Validation_CTh=table2array(Field_CTh);
Validation_SedTrap=table2array(Field_SedTrap);

%Optional plotting
if toplot==1
    
    figure('Position',[50 100 800 600])
    subplot(4,1,1)
    pcolor(time,z,Phy_track)
    shading flat
    colorbar
    axis ij
    caxis([0 40])
    colormap(turbo)
    title('Model Phy')
    
    subplot(4,1,2)
    pcolor(time,z,Phy_hold)
    shading flat
    colorbar
    axis ij
    caxis([0 40])
    colormap(turbo)
    title('Measured Phy')
    
    subplot(4,1,3)
    pcolor(time,z,Phy_track-Phy_hold)
    shading flat
    colorbar
    axis ij
    caxis([-0.01 0.01])
    title('Difference')
    
    subplot(4,1,4)
    pcolor(time,z,(Phy_track-Phy_hold)./Phy_hold*100)
    shading flat
    colorbar
    axis ij
    caxis([-3 3])
    title('Difference (Percent Error)')
    

    figure('Position',[50 50 800 600])
    subplot(3,1,1)
    pcolor(time,z,Det_track)
    shading flat
    colorbar
    axis ij
    caxis([0 40])
    colormap(turbo)
    title('Model Det')
    
    subplot(3,1,2)
    pcolor(time,z,Det_hold)
    shading flat
    colorbar
    axis ij
    caxis([0 40])
    colormap(turbo)
    title('Measured Det')
    
    subplot(3,1,3)
    pcolor(time,z,Det_track-Det_hold)
    shading flat
    colorbar
    axis ij
    caxis([-0.01 0.01])
    title('Difference')

    
    figure('Position',[350 50 800 600])
    ax(1)=subplot(3,1,1);
    pcolor(time,z,dTh_track+PTh_track)
    shading flat
    colorbar
    axis ij
    title('Total Thorium (dpm L^-^1)')
    
    ax(2)=subplot(3,1,2);
    pcolor(time,z,PTh_track./(dTh_track+PTh_track))
    shading flat
    colorbar
    axis ij
    title('Fraction Particulate Thorium')
    
    
    ax(3)=subplot(3,1,3);
    pcolor(time,z,POC_track./PTh_track)
    shading flat
    colorbar
    axis ij
    title('Carbon:Thorium (umol C dpm^-^1)')
    hold on
    plot(Validation_CTh(:,1),Validation_CTh(:,2),'ok','MarkerFaceColor','k','MarkerSize',8)
    scatter(Validation_CTh(:,1),Validation_CTh(:,2),40,Validation_CTh(:,5),'filled')
    colormap(ax(3),turbo)
    
    
    figure('Position',[550 50 800 300])
    plot(time,SumSink_track./SumThSink_track,'-b','Color',[0 0.5 1],'LineWidth',2)
    hold on
    for i=1:length(Validation_SedTrap(:,1))
        plot([Validation_SedTrap(i,3) Validation_SedTrap(i,4)],[Validation_SedTrap(i,11) Validation_SedTrap(i,11)],'-r','Color',[1 0 0],'LineWidth',2)
    end
    ylabel('C:Th Ratio (umol dpm^-^1)')
    title('Sinking Particles')
end

for i=1:length(Validation_CTh(:,1))
    Validation_CTh(i,7)=POC_track(Validation_CTh(i,4),Validation_CTh(i,3))/PTh_track(Validation_CTh(i,4),Validation_CTh(i,3));
end
norm_misfit_CTh = (Validation_CTh(:,7)-Validation_CTh(:,5))./Validation_CTh(:,6);  %Error divided by uncertainty

for i=1:length(Validation_SedTrap(:,1))
    %Total sediment trap
    tempPOC = mean(SumSink_track(Validation_SedTrap(i,1):Validation_SedTrap(i,2)));
    tempTh = mean(SumThSink_track(Validation_SedTrap(i,1):Validation_SedTrap(i,2)));
    Validation_SedTrap(i,19)=tempPOC/tempTh;
    
    %<200
    tempPOC = mean(POCSink_track(Validation_SedTrap(i,1):Validation_SedTrap(i,2)));
    tempTh = mean(ThSink_track(Validation_SedTrap(i,1):Validation_SedTrap(i,2)));
    Validation_SedTrap(i,20)=tempPOC/tempTh;
    
    %>200
    tempPOC = mean(FPCFlux_track(Validation_SedTrap(i,1):Validation_SedTrap(i,2)));
    tempTh = mean(FPThFlux_track(Validation_SedTrap(i,1):Validation_SedTrap(i,2)));
    Validation_SedTrap(i,21)=tempPOC/tempTh;
    
    
end
norm_misfit_SedTrap = (Validation_SedTrap(:,11)-Validation_SedTrap(:,19))./Validation_SedTrap(:,12);
norm_misfit_SedTrapSmall = (Validation_SedTrap(:,7)-Validation_SedTrap(:,20))./Validation_SedTrap(:,8);
norm_misfit_SedTrapLarge = (Validation_SedTrap(:,9)-Validation_SedTrap(:,21))./Validation_SedTrap(:,10);

norm_misfit = [norm_misfit_CTh;norm_misfit_SedTrap;norm_misfit_SedTrapSmall;norm_misfit_SedTrapLarge];

norm_misfit(find(isnan(norm_misfit)))=[];