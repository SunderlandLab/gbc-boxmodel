%==========================================================================
% OBJECTIVE
%   Process historical anthropogeic Hg emissions.
%
% REFERENCES
%   Amos, H. M., et al. (2013), Legacy impacts of all-time anthropogenic 
%     emissions on the global mercury cycle, Glob. Biogeochem. Cycle, 
%     27(2), 410-421.
%   Engstrom, D. R., et al. (2014), Atmospheric hg emissions from 
%     preindustrial gold and silver extraction in the americas: A 
%     reevaluation from lake-sediment archives, Environ. Sci. Technol., 
%     48(12), 6533-6543.
%   Horowitz, H. M., et al. (2014), Historical mercury releases from 
%     commercial products: Global environmental implications, Environ. Sci.
%     Technol., 48(17), 10242-10250.
%   Nriagu, J. O. (1993), Legacy of mercury pollution, Nature, 363(6430), 
%     589-589.
%   Streets, D. G., et al. (2009), Projections of global mercury emissions 
%     in 2050, Environ. Sci. Technol., 43(8), 2983-2988.
%   Streets, D. G., et al. (2011), All-time releases of mercury to the 
%     atmosphere from human activities, Environ. Sci. Technol., 45(24), 
%     10485-10491.   
%
% REVISION HISTORY
%   02 Nov 2011 - HMA - now interpolate with pchip() instead of spline()
%   08 Nov 2011 - HMA - use built-in MATLAB functions to optimize the code
%   14 Jan 2012 - hma - extrend pre-1850 anthropogenic emissions to 2000 BC
%   17 Jan 2012 - hma - add time-dependent volcanic eruptions
%   27 Apr 2012 - hma - add David Streets' new 2050 SRES scenarios
%   03 May 2012 - hma - extrapolate future scenarios to 2100
%   30 Jul 2012 - hma - clean up code/comments
%   30 Jan 2014 - HMA - update filepath to Streets anthro emission files
%   02 May 2014 - HMA - create figure of primary emissions from 1500 to
%                       2050 for private defense slides
%   02 Sep 2014 - HMA - add option to decrease historical mining emissions
%                       (Engstrom et al., 2014)
%   08 Sep 2014 - HMA - clean up code and comments for public release
%==========================================================================

%%
%--------------------------------------------------------------------------
% Emission inventory from Streets et al. (2011)
%--------------------------------------------------------------------------

% read in emissions from .txt file
load('AnthroEmissAllTime_20120112.txt')

% parse data
Syear   = AnthroEmissAllTime_20120112(:,1); % 2000 BC to 2008 AD, decadal
Streets = AnthroEmissAllTime_20120112(:,2); % Mg/yr

Streets2 = AnthroEmissAllTime_20120112(:,4); % 80% upper confidence interval for emissions
Streets3 = AnthroEmissAllTime_20120112(:,3); % 80% lower confidence interval for emissions


%%
%--------------------------------------------------------------------------
% Option to modify Streets et al. (2011) inventory by reducing historical
% mining sources, as suggested by Engstrom et al. (2014)
%--------------------------------------------------------------------------

% By the end of the 1920s cyanide had completely supplanted Hg for Au and 
% Ag extraction, so apply 50% reduction to all decades before the 1920s
if strcmp(emissInven,'Engstrom');
    % global totals for with mining sources decreased by 50%
    Streets( Syear < 1850  ) = 0.5*Streets( Syear < 1850  );
    Streets( Syear == 1850 ) = Streets( Syear == 1850 ) - 495/2; % where 495 Mg/yr is the total of Hg, Au, and Ag mining
    Streets( Syear == 1860 ) = Streets( Syear == 1860 ) - 620/2;
    Streets( Syear == 1870 ) = Streets( Syear == 1870 ) - 1491/2;
    Streets( Syear == 1880 ) = Streets( Syear == 1880 ) - 2157/2;
    Streets( Syear == 1890 ) = Streets( Syear == 1890 ) - 2489/2;
    Streets( Syear == 1900 ) = Streets( Syear == 1900 ) - 2140/2;
    Streets( Syear == 1910 ) = Streets( Syear == 1910 ) - 1472/2;
end

%%

% interpolate to annual resolution
Time    = (Syear(1):Syear(end));
Anthro  = pchip(Syear, Streets , Time);  % mid
Anthro2 = pchip(Syear, Streets2, Time);  % high
Anthro3 = pchip(Syear, Streets3, Time);  % low

% best estimate
basecase = Anthro; % dummy


%%
%--------------------------------------------------------------------------
% Future emission projections from Amos et al. 2013, based in part on
% Streets et al., 2009
%--------------------------------------------------------------------------

% David Streets' revised 2050 SRES scenarios
clear n; % for safety's sake
FTime = [Time, 2009:1:2050]; % future time, annual resoltuion

% initialize future scenarios with 2009
% AnthroA1B  = [Anthro, Anthro(end) + 55.91];   % A1B in 2009
AnthroA1B  = [basecase, Anthro(end) + 55.91];   % A1B in 2009

% loop thru 2010 - 2015
clear n; % for safety's sake
n = length(AnthroA1B) + 1;
for j = n:n+5; % business as usual until 2015
%for j = n:n+9; % business as usual until 2019
%for j = n:n+19;% business as usual until 2029    
    AnthroA1B(j) = AnthroA1B(j-1) + 55.91;   % A1B
end

% all scenarios are business-as-usual between 2008-2015 and then in 2016 we
% assume that the UNEP treaty is implemented
AnthroB1   = AnthroA1B;   % B1
AnthroBest = AnthroA1B;   % best case
AnthroZero = AnthroA1B;   % zero emissions

% implement scenarios beginning in 2016
clear j n m;
n = length(AnthroA1B)+1;
m = 2050 - 2016; % treaty comes into force in 2016
%m = 2050 - 2020; % treaty comes into force in 2020
%m = 2050 - 2030; % treaty comes into force in 2030
for j = n:n+m
    AnthroA1B(j)  = AnthroA1B(j-1)  + 55.91; % A1B
    AnthroB1(j)   = AnthroB1(j-1)   + 0    ; % constant
    AnthroBest(j) = AnthroBest(j-1) - 38.6;  % best case scenario
    AnthroZero(j) = 0;                       % zero emissions as of 2016    
end

%%
%--------------------------------------------------------------------------
% Emissions from intentional use of Hg in commericial products and
% processes
%
% Reference: Horowitz et al. (2014)
%--------------------------------------------------------------------------

% If you use the Horowitz et al. (2014) inventory, this automatically sets
% a logical that tells the code to add commericial Hg emissions to the
% Streets et al. (2011) inventory
if strcmp(emissInven,'Horowitz');
    Lprod = 1;
    forWeb_product_emissions
else
    Lprod = 0;
end


%--------------------------------------------------------------------------
% PLOTS
%--------------------------------------------------------------------------

if Lplot;

    % primary anthropogenic emissions, including future scenarios
    hfig = figure(70);
    set(hfig,'units','normalized','Position',[0.1 0.4 0.5 .7])
    set(gcf,'Color',[1 1 1])
    subplot(2,1,1)
    set(gca,'FontSize',18)
    hold on;
    plot(FTime,AnthroB1  ,'b','LineWidth',3)
    plot(FTime,AnthroBest,'g','LineWidth',3)
    plot(FTime,AnthroZero,'c','LineWidth',3)
    plot(FTime,AnthroA1B ,'r','LineWidth',3)
    plot(Time ,basecase  ,'k','LineWidth',4)
    xlabel('Year (AD)')
    ylabel('(Mg a^{-1}) ')
    xlim([1500 2050])
    title('Primary Anthropogenic Emissions (Streets 2011)')
    hold off;

     % cummulative anthropogenic emissions
    subplot(2,1,2)
    set(gca,'FontSize',18)
    hold on;
    plot(FTime,1e-3*cumsum(AnthroB1)  ,'b','LineWidth',3)
    plot(FTime,1e-3*cumsum(AnthroBest),'g','LineWidth',3)
    plot(FTime,1e-3*cumsum(AnthroZero),'c','LineWidth',3)
    plot(FTime,1e-3*cumsum(AnthroA1B) ,'r','LineWidth',3)
    plot(Time ,1e-3*cumsum(basecase)    ,'k','LineWidth',4)
    xlabel('Year (AD)')
    ylim([0 5e2])
    ylabel('(Gg) ')
    legend('Constant','Mercury controls','Zero','A1B','Historical','Location','NorthWest')
    xlim([1500 2050])
    title('Cumulative Emissions')
    hold off;
    
end
