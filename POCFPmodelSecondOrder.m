function [norm_misfit] = POCFPmodelSecondOrder(k_sorp,k_desorp,Eg_Th,starttimes,endtimes,toplot)

%Eg_Th = Fraction of thorium consumed by krill that passes to their fecal pellets

% clearvars
% close all
% k_sorp = 1/30;
% k_desorp = 1/365;
% Eg_Th = 0.5;
% starttimes = [9.5,50];
% endtimes = [40,137];
% toplot = 1;

Eg_C = 0.3;  %Fraction of carbon consumed by prill that passes to their fecal pellets (i.e., egestion efficiency = 1 - assimilation efficiency)

%Diagnostics
CompileValidationData=0;

%Setting up boundaries and forcing
load('POC-FP model inputs 1.mat')
load('POC-FP model inputs 2.mat')
load('POC-FP model inputs 3.mat')
POC_hold = POC(1:find(z==65),:);
PTh_hold = PTh(1:find(z==65),:);
dTh_hold = Th(1:find(z==65),:)-PTh_hold;
dt = time(2)-time(1);
C14PP = C14PP(1:find(z==65),:);
Grazing = Grazing(1:find(z==65),:);
Remin = Remin(1:find(z==65),:);
omega = omega(1:find(z==65),:);
Kz = Kz(1:find(z_edge==67.5),:);
z=z(1:find(z==65));
z_edge=z_edge(1:find(z_edge==67.5));
z_thick=z_edge(2:end)-z_edge(1:end-1);

%Unchanging parameters
lam = 0.0287612937991679;


%Initializing variables
POC_track = NaN(size(POC_hold));
PTh_track = NaN(size(POC_hold));
dTh_track = NaN(size(POC_hold));
POCSink_track = NaN(size(time));
ThSink_track = NaN(size(time));
FPCFlux_track = NaN(size(time));
FPThFlux_track = NaN(size(time));

%Model Run - First portion
POC = POC_hold(:,find(time==starttimes(1)));
PTh = PTh_hold(:,find(time==starttimes(1)));
dTh = dTh_hold(:,find(time==starttimes(1)));
clear Th_temp
for i=find(time==starttimes(1)):find(time==endtimes(1))

    %Sinking
    [POC_out,POCFlux] = sinking(POC,z_edge',omega(:,i),dt);
    POC = POC_out;
    
    [PTh_out,PThFlux] = sinking(PTh,z_edge',omega(:,i),dt);
    PTh = PTh_out;
    
    %Mixing
    [Coeff0,Coeff1,Coeff2,BottomCoeff]=CalculateMixingCoefficients(POC,z,z_edge,Kz(:,i),dt,POC_hold(end,i));
    [POC_out]=mixing_ftcs(POC,POC_hold(end,i),Coeff0,Coeff1,Coeff2);
    POC = POC_out;
    
    [PTh_out]=mixing_ftcs(PTh,PTh_hold(end,i),Coeff0,Coeff1,Coeff2);
    PTh = PTh_out;
    
    [dTh_out]=mixing_ftcs(dTh,dTh_hold(end,i),Coeff0,Coeff1,Coeff2);
    dTh = dTh_out;

    %Primary Production
    POC = POC + C14PP(:,i)*dt;
    
    %Grazing
    PTh = PTh - Grazing(:,i)*dt .* PTh./POC;
    dTh = dTh + Grazing(:,i)*dt*(1-Eg_Th) .* PTh./POC;
    POC = POC - Grazing(:,i)*dt;
    FPCFlux_track(i) = sum(Grazing(:,i)*dt.*z_thick')*Eg_C;
    FPThFlux_track(i) = sum(Grazing(:,i)*dt.* PTh./POC.*z_thick')*Eg_Th;
    
    %Remineralization
    PTh = PTh - Remin(:,i)*dt .* PTh./POC;
    dTh = dTh + Remin(:,i)*dt .* PTh./POC;
    POC = POC - Remin(:,i)*dt;
    
    %Sorption
    PTh = PTh + k_sorp*POC.*dTh*dt;
    dTh = dTh - k_sorp*POC.*dTh*dt;
    
    %Desorption
    PTh = PTh - k_desorp*PTh*dt;
    dTh = dTh + k_desorp*PTh*dt;
    
    %Th-234 Decay and production
    PTh = PTh - PTh*lam*dt;
    dTh = dTh - dTh*lam*dt;
    dTh = dTh + U238(:,i)*lam*dt;
    
    POC_track(:,i+1)=POC;
    PTh_track(:,i+1)=PTh;
    dTh_track(:,i+1)=dTh;
    ThSink_track(i)=PThFlux(find(z==50));
    POCSink_track(i)=POCFlux(find(z==50));
end

%Model Run - Second portion
POC = POC_hold(:,find(time==starttimes(2)));
PTh = PTh_hold(:,find(time==starttimes(2)));
Th_temp = Th(1:find(z==65),find(time==starttimes(2)));
dTh = Th_temp - PTh;
clear Th_temp
for i=find(time==starttimes(2)):find(time==endtimes(2))

    %Sinking
    [POC_out,POCFlux] = sinking(POC,z_edge',omega(:,i),dt);
    POC = POC_out;
    
    [PTh_out,PThFlux] = sinking(PTh,z_edge',omega(:,i),dt);
    PTh = PTh_out;
    
    %Mixing
    [Coeff0,Coeff1,Coeff2,BottomCoeff]=CalculateMixingCoefficients(POC,z,z_edge,Kz(:,i),dt,POC_hold(end,i));
    [POC_out]=mixing_ftcs(POC,POC_hold(end,i),Coeff0,Coeff1,Coeff2);
    POC = POC_out;
    
    [PTh_out]=mixing_ftcs(PTh,PTh_hold(end,i),Coeff0,Coeff1,Coeff2);
    PTh = PTh_out;
    
    [dTh_out]=mixing_ftcs(dTh,dTh_hold(end,i),Coeff0,Coeff1,Coeff2);
    dTh = dTh_out;

    %Primary Production
    POC = POC + C14PP(:,i)*dt;
    
    %Grazing
    PTh = PTh - Grazing(:,i)*dt .* PTh./POC;
    dTh = dTh + Grazing(:,i)*dt*(1-Eg_Th) .* PTh./POC;
    POC = POC - Grazing(:,i)*dt;
    FPCFlux_track(i) = sum(Grazing(:,i)*dt.*z_thick')*Eg_C;
    FPThFlux_track(i) = sum(Grazing(:,i)*dt.* PTh./POC.*z_thick')*Eg_Th;
    
    %Remineralization
    PTh = PTh - Remin(:,i)*dt .* PTh./POC;
    dTh = dTh + Remin(:,i)*dt .* PTh./POC;
    POC = POC - Remin(:,i)*dt;
    
    %Sorption
    PTh = PTh + k_sorp*POC.*dTh*dt;
    dTh = dTh - k_sorp*POC.*dTh*dt;
    
    %Desorption
    PTh = PTh - k_desorp*PTh*dt;
    dTh = dTh + k_desorp*PTh*dt;
    
    %Th-234 Decay and production
    PTh = PTh - PTh*lam*dt;
    dTh = dTh - dTh*lam*dt;
    dTh = dTh + U238(:,i)*lam*dt;
    
    POC_track(:,i+1)=POC;
    PTh_track(:,i+1)=PTh;
    dTh_track(:,i+1)=dTh;
    ThSink_track(i)=PThFlux(find(z==50));
    POCSink_track(i)=POCFlux(find(z==50));
end

%Combining sinking POC and fecal pellets
SumSink_track = POCSink_track + FPCFlux_track;
SumThSink_track = ThSink_track + FPThFlux_track;

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
    
    
    figure('Position',[350 50 800 900])
    ax(1)=subplot(4,1,1);
    pcolor(time,z,dTh_track+PTh_track)
    shading flat
    colorbar
    axis ij
    title(['POC-FP-Model 2^n^d Order',char(10),'Total Thorium (dpm L^-^1)'])
    
    ax(2)=subplot(4,1,2);
    pcolor(time,z,PTh_track./(dTh_track+PTh_track))
    shading flat
    colorbar
    axis ij
    title('Fraction Particulate Thorium')
    
    
    ax(3)=subplot(4,1,3);
    pcolor(time,z,POC_track./PTh_track)
    shading flat
    colorbar
    axis ij
    title('Carbon:Thorium (umol C dpm^-^1)')
    hold on
    plot(Validation_CTh(:,1),Validation_CTh(:,2),'ok','MarkerFaceColor','k','MarkerSize',8)
    scatter(Validation_CTh(:,1),Validation_CTh(:,2),40,Validation_CTh(:,5),'filled')
    colormap(ax(3),turbo)
    
    
    ax(4)=subplot(4,1,4);
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