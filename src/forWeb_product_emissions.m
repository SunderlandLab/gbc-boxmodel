%==========================================================================
% OBJECTIVE
%   Inputs of mercury (Hg) to landfills, terrestrial reservoirs, and
%   atmosphere from intentional uses 1850 - 2008
%   
% REFERENCE
%   Horowitz, H. M., et al. (2014), Historical mercury releases from 
%     commercial products: Global environmental implications, Environ. Sci.
%     Technol., 48(17), 10242-10250.
% 
% REVISION HISTORY
%   11 Jun 2013 - HMH - add different emissions option from intentional
%                       uses after reducing global Hg consumption in paint
%   28 Dec 2013 - HMA - inherited code from HMH. Clean up to merge with EMS's
%                       updated ocean code. 
%                       Move all flags/logicals to main.m.
%   18 Mar 2014 - HMA - updated emissions from Hannah 
%   08 Sep 2014 - HMA - clean up code and comments for public release
%
% Hannah M. Horowitz
% hmhorow@fas.harvard.edu ; hmhorow@gmail.com ; hmhorow@post.harvard.edu
%==========================================================================

% Product releases from Horowitz et al., 2014
% skipping 1st row which is column headers
prod_emissions = csvread('product_emissions_031114_corrected_to_2008_only.csv',1,0); % HMA, 18 Mar 2014
disp('Reading file: product_emissions_031114_corrected_to_2008_only.csv')

% parse output
P_atm           = prod_emissions(:,1); % this is only the additional atmospheric releases from my inventory 
P_w             = prod_emissions(:,2);
P_soil          = prod_emissions(:,3);
P_lf            = prod_emissions(:,4);
E_waste         = prod_emissions(:,5); % needs to be subtracted from Streets'atmospheric emissions inventory to avoid double counting with my additional releases
Streets_overlap = prod_emissions(:,6) - prod_emissions(:,1); % emissions from chlor alkali, ASGM, and large-scale gold and silver mining

% since we need to go from 2000 BC to 2008 AD, we will need to add
% additional elements at the beginning to be all 0s before 1850

n_el      = length(Syear);  % total # of elements for 2000 BC - 2008 AD
n_prod    = length(P_atm);  % # of elements for just 1850 - 2008 AD
extra     = n_el - n_prod;  % the number of elements that must be added to the product Hg for 2000 BC - 1850 AD

pre_1850  = zeros(extra,1); % column of 0s from 2000 BC through 1840

% add to each column of 1850 - 2008 data
P_atm           = [pre_1850; P_atm          ]; 
P_w             = [pre_1850; P_w            ]; 
P_soil          = [pre_1850; P_soil         ]; 
P_lf            = [pre_1850; P_lf           ]; 
E_waste         = [pre_1850; E_waste        ]; 
Streets_overlap = [pre_1850; Streets_overlap]; 

% interpolate each to annual resolution
Time              = (Syear(1):Syear(end));  % time increasing by 1 year instead of every 10 years
Ep_atm            = pchip(Syear, P_atm           , Time); 
Ep_w              = pchip(Syear, P_w             , Time); 
Ep_soil           = pchip(Syear, P_soil          , Time); 
Ep_lf             = pchip(Syear, P_lf            , Time); 
E_Streets_waste   = pchip(Syear, E_waste         , Time); 
E_Streets_overlap = pchip(Syear, Streets_overlap , Time); 


% plot all the emissions for each reservoir
% manual stacking of lines  by summing each one to the previous ones
if Lplot;
    figure(40)
    set(gca, 'FontSize',18)
    hold on
    h1=plot(Time,Ep_atm,'r','LineWidth',3);
    h2=plot(Time,Ep_w+Ep_atm,'b','LineWidth',3);
    h3=plot(Time,Ep_soil+Ep_w+Ep_atm,'g','LineWidth',3);
    h4=plot(Time,Ep_lf+Ep_soil+Ep_w+Ep_atm,'k','LineWidth',3);
    legend([h1 h2 h3 h4],'Atmosphere','Water','Soil','Landfills','Location','Best');
    xlim([1850 2008])
    title ('Additional releases from intentional uses of Hg')
    set(gca, 'FontSize',18)
    xlabel('Time (years)')
    ylabel('Mg of Hg')
    hold off
end

% create area plot
all_emissions = [Ep_atm' Ep_w' Ep_soil' Ep_lf'];

% plot all the emissions for each reservoir
if Lplot;
    figure(41)
    set(gca, 'FontSize',18)
    area(Time, all_emissions)
    colormap jet
    xlabel('Year (AD)')
    ylabel('Hg (Mg a^{-1})')
    title ('Additional releases from intentional uses of Hg')
    ylim([0 7000])
    xlim([1850 2008])
    legend('Atmosphere','Water','Soil','Landfills','Location','Best')
end

