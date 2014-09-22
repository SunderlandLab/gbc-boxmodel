%==========================================================================
% OBJECTIVE
%   This script calculates rate constants as first order linear processes
%   from the present global Hg budget.  We assume that rate constants are 
%   the same between the preindustrial and present-day period and that 
%   biomass burning is neglible for the preindustrial simulation (Selin et 
%   al. 2008)
%
% REFERENCES
%    Amos, H. M., et al. (2014), Global biogeochemical implications of 
%      mercury discharges from rivers and sediment burial, Environ. Sci.
%      Technol., 48(16), 9514-9522.
%    Bagnato, E., et al. (2011), New clues on the contribution of earth's 
%      volcanism to the global mercury cycle, Bull. Volcanol., 73(5), 497-
%      510.
%    Hararuk, O., et al. (2013), Modeling the sensitivity of soil mercury 
%      storage to climate-induced changes in soil carbon pools, 
%      Biogeosciences, 10, 2393?2407.
%    Holmes, C. D., et al. (2010), Global atmospheric model for mercury 
%      including oxidation by bromine atoms, Atmos. Chem. Phys., 10, 12037
%      -12057.
%    Pirrone, N., et al. (2010), Global mercury emissions to the atmosphere
%      from anthropogenic and natural sources, Atmos. Chem. Phys., 10(13), 
%      5951-5964.
%    Selin, N. E., et al. (2008), Global 3-d land-ocean-atmosphere model 
%      for mercury: Present-day versus preindustrial cycles and 
%      anthropogenic enrichment factors for deposition Global Biogeochem. 
%      Cycles, 22(3), GB3099.
%    Smith-Downey, N. V., et al. (2010), Anthropogenic impacts on global 
%      storage and emissions of mercury from terrestrial soils: Insights 
%      from a new global model, J. Geophys. Res.-Biogeosci., 115, G03008.
%   Soerensen, A. L., et al. (2010), An improved global model for air-sea 
%      exchange of mercury: High concentrations over the north atlantic, 
%      Environ. Sci. Technol., 44(22), 8574-8580.   
%
% MODIFCATION HISTORY
%   05 Jun 2011 - EMS - last dated modification
%   08 Nov 2011 - HMA - clean up comments. Add variables for fluxes and
%                       explicity calculate rate constants, so that the k
%                       values aren't just hard coded. 
%   09 Nov 2011 - HMA - correct the gross ocean fluxes (using Anne's
%                       diagnostics fix for GEOS-Chem)
%   06 Dec 2011 - hma - specify different biomass burning emissions for the
%                       pre-1450 period and 1450-2008
%   07 Dec 2011 - hma - fix error in biomass burning, should be constant
%                       not f(t)
%   19 Jul 2012 - hma - make geogenic emissions and external forcing
%   30 Jul 2012 - hma - clean up code/comments
%   19 Mar 2014 - hma - Correct gross Hg0 uptake and evasion.
%                       Correct Dep_oHgII from 3900 to 3600 Mg/yr based on
%                       Holmes et al. 2010.
%   10 Apr 2014 - hma - add LBessTerr logical to test effect of using Bess'
%                       updated numbers for terrestrial cycling
%   08 Sep 2014 - HMA - clean up code and comments for public release
%==========================================================================


%--------------------------------------------------------------------------
% Reservoirs estimates for present day (Mg)
%--------------------------------------------------------------------------

% Atmosphere
Ratm       = 5000;               % Holmes et al.(2010)

% Ocean
Rocs       = 2910;               % Soerensen et al. (2010)
Roci       = 134000;             % Sunderland and Mason (2007)
Rocd       = 220649;             % Sunderland and Mason (2007)

%  Terrestrial
Rtf        = 9620;               % Leaf, fast and intermediate pools from Smith-Downey et al (2010)
Rts        = 34900;              % Slow pool from Smith-Downey et al. (2010)
Rta        = 193600;             % Armored pool from Smith-Downey et al. (2010)

% Deep Mineral Reserves
Rm         = 3e11;               % Deep reservoirs (Selin et al. 2008)

%--------------------------------------------------------------------------
% Fluxes for present day (Mg/year)
%--------------------------------------------------------------------------

% Atmosphere
Dep_oHgII  = 3600;%3900;      % Hg(II) deposition to the ocean (Holmes et al., 2010), corrected from 3900 to 3600 Mg/yr by HMA 19 Mar 2014
Dep_tHgII  = 1500;            % Hg(II) deposition to land (Holmes et al., 2010)
Dep_tHg0   = 1500;            % Hg(0) deposition to land  (Holmes et al., 2010)
                              
% Hg0 air - sea exchange
netEv_Hg0  = 3000;            % net evasion from surface ocean to atmosphere (Soerensen et al., 2010)
Ev_Hg0_ocs = 4700;            % gross ocean evasion to the atmosphere (Soerensen et al., 2010)
Upt_oHg0   = Ev_Hg0_ocs - netEv_Hg0; % gross uptake of Hg0 from the atmopshere (1700 Mg/yr from Holmes et al., 2010)

% Surface ocean
ps_ocs     = 3320;            % particle settling 
vert_ocsi  = 5100;            % gross detrainment flux, surface to intermediate (ref: HMA, v9-01-02, 2005 avg)

% Intermediate ocean
ps_oci     = 480;             % particle settling             
vert_ocis  = 7100;            % vertical seawater flow, intermediate to surface
vert_ocid  = 335;             % vertical seawater flow, intermediate to deep

% Deep ocean
ps_ocd     = 210;             % particle settling, burial in deep sediments
vert_ocdi  = 175;             % vertical sea water flow, deep to intermediate

% biomass burning
E_bb_1     = 0;               % pre-anthropogenic biomass burning emissions 
E_bb_2     = 300;             % anthropogenic biomass burning emissions (Holmes et al., 2010)

% rivers and biomass burning: assuming 75% vegetation (all fast) + 25%
% soils (fast,slow,armored) partitioned by C storage
fveg       = 0.95;            % fraction to vegetation
fsoil      = 0.05;            % fraction to soil

% fraction of carbon in each pool (Smith-Downey et al., 2010)
fCfast     = 0.2185;          % fast reservoir
fCslow     = 0.5057;          % slow reservoir
fCarmored  = 0.2759;          % armored reservoir

% Fast terrestrial
E_Te_rf    = 460;             % evasion due to respiration of organic carbon
E_Te_p     = 845;             % photoreduction and re-release of deposited Hg0 (Smith-Downey et al. 2010)
Te_exfs    = 325;             % exchange among soil pools, fast pool to slow pool
Te_exfa    = 9;               % exchange among soil pools, fast pool to armored pool

% Slow terrestrial
E_Te_rs    = 250;             % evasion due to respiration of organic carbon
Te_exsf    = 205;             % exchange among soil pools, slow pool to fast pool
Te_exsa    = 0.5;             % exchange among soil pools, slow pool to armored pool

% Armored terrestrial
E_Te_ra    = 25;              % evasion due to respiration of organic carbon
Te_exaf    = 15;              % exchange among soil pools, armored pool to fast pool
Te_exam    = 0;               % exchange from armored pool to mineral pool 

% Deep mineral reservoir 
E_geo      = 90;              % geogenic emissions (Pironne et al., 2010; Bagnato et al., 2011)

%--------------------------------------------------------------------------
% Atmospheric rates (1/year)
%--------------------------------------------------------------------------
k_A_oHgII  = Dep_oHgII / Ratm;   % HgII deposition to surface ocean; min = 0.7200 (CDH); max = 0.8913 (ESC)
k_A_tHgII  = Dep_tHgII / Ratm;   % HgII deposition to terrestrial surfaces;  min = 0.2927 (NES); max = 0.3043 (ESC)
k_A_tHg0   = Dep_tHg0  / Ratm;   % Hg0  deposition to terrestrial surfaces; 
                                 %   Assumed to represent terrestrial leaf uptake; 
                                 %   min = 0.2927 (NES); max = 0.3478 (ESC)
k_A_oHg0   = Upt_oHg0  / Ratm;   % gross Hg0 uptake by the surface ocean

% fraction of atmopsheric deposition to...
fdep_tf    = 0.5027;             % the fast soil pool
fdep_ts    = 0.3213;             % the slow soil pool
fdep_ta    = 0.1753;             % the fast armored pool
  

%--------------------------------------------------------------------------
% Surface ocean rates (1/year)
%--------------------------------------------------------------------------    
k_Oc_ev    = Ev_Hg0_ocs / Rocs;  % evasion Hg0 (gross flux); min = 0.7400 (CDH); max = 1.8485 (ESC) 
k_Oc_sp1   = ps_ocs     / Rocs;  % particle settling 
k_Oc_vsi   = vert_ocsi  / Rocs;  % gross detrainment to intermediate ocean

%--------------------------------------------------------------------------
% Intermediate ocean rates (1/year) 
%--------------------------------------------------------------------------
k_Oc_sp2   = ps_oci    / Roci;   % particle settling    
k_Oc_vis   = vert_ocis / Roci;   % vertical seawater flow, intermediate to surface
k_Oc_vid   = vert_ocid / Roci;   % vertical seawater flow, intermediate to deep

%--------------------------------------------------------------------------
% Deep ocean rates (1/year)
%--------------------------------------------------------------------------
k_Oc_sp3   = ps_ocd    / Rocd;   % particle settling, deep ocean burial in deep sediments 
k_Oc_vdi   = vert_ocdi / Rocd;   % vertical seawater flow, deep to intermediate

%--------------------------------------------------------------------------
% Fast terrestrial reservoir rates (1/year)
%--------------------------------------------------------------------------

% Includes vegetation, ice, and fast+intermediate carbon pools from 
% Smith-Downey et al. (2010)

k_Te_rf    = E_Te_rf / Rtf;      % respiration
k_Te_p     = E_Te_p  / Rtf;      % photoreduction and re-release of deposited Hg0
k_T_exfs   = Te_exfs / Rtf;      % exchange among soil pools, fast pool to slow pool
k_T_exfa   = Te_exfa / Rtf;      % exchange among soil pools, fast pool to armored pool 

% biomass burning 
k_Te_BBf_1 = (E_bb_1*fveg   + (E_bb_1*fsoil*fCfast))   / Rtf;  % pre-1450
k_Te_BBf_2 = (E_bb_2*fveg   + (E_bb_2*fsoil*fCfast))   / Rtf;  % 1450-2000

%--------------------------------------------------------------------------
% Slow terrestrial reservoir rates (1/year)
%--------------------------------------------------------------------------

k_Te_rs    = E_Te_rs / Rts;       % evasion due to respiration of organic carbon
k_T_exsf   = Te_exsf / Rts;       % exchange among soil pools, slow pool to fast pool
k_T_exsa   = Te_exsa / Rts;       % exchange among soil pools, slow pool to armored pool

% biomass burning 
k_Te_BBs_1 = (E_bb_1*fsoil*fCslow)   / Rts;   % pre-1450
k_Te_BBs_2 = (E_bb_2*fsoil*fCslow)   / Rts;   % 1450-2008

%--------------------------------------------------------------------------
% Armored terrestrial reservoir rates (1/year)
%--------------------------------------------------------------------------

k_Te_ra    = E_Te_ra / Rta;       % evasion due to respiration of organic carbon
k_T_exaf   = Te_exaf / Rta;       % exchange among soil pools, armored pool to fast pool
k_T_exam   = Te_exam / Rta;       % exchange from armored pool to mineral pool (18 Dec 2011, hma)

% biomass burning 
k_Te_BBa_1   = (E_bb_1*fsoil*fCarmored)   / Rta; % pre-1450
k_Te_BBa_2   = (E_bb_2*fsoil*fCarmored)   / Rta; % 1450-2008


% Decrease soil Smith-Downey 2010 re-emissions from respiration and
% photoreduction to be more in line with new observations and estimates
% from Hararuk, Obris et al. (2013) and USGS. The adjustment factor may
% need to be revisted as more data is published. 
reduce_factor   = 1/10;  

% fast terrestrial
k_Te_rf = reduce_factor*k_Te_rf;
k_Te_p  = reduce_factor*k_Te_p;

% slow soil
k_Te_rs = reduce_factor*k_Te_rs;

% armored soil
k_Te_ra = k_Te_ra*reduce_factor;


%--------------------------------------------------------------------------
% Rivers
%--------------------------------------------------------------------------

% total discharged to ocean margins
Te_riv_margin  = IHgD_pristine + IHgP_pristine;

% global fraction of riverine HgP reaching the open oceans (Table 2, Amos et al. 2014)
if strcmp(Lriver_FHgP,'Walsh');
    f_HgPexport    = 0.30;
elseif strcmp(Lriver_FHgP,'Chester');
    f_HgPexport    = 0.10;
end

% total reaching the open ocean
Te_riv_ocean   = IHgD_pristine + f_HgPexport*IHgP_pristine;

% First-order rate coefficients (1/yr)

% Riverine discharge of terrestrial Hg to ocean margins
k_T_riv_f      = (Te_riv_margin*fveg + (Te_riv_margin*fsoil*fCfast)) / Rtf;  % fast
k_T_riv_s      = (Te_riv_margin*fsoil*fCslow) / Rts;                         % slow
k_T_riv_a      = (Te_riv_margin*fsoil*fCarmored) / Rta;                      % armored

% Of the riverine discharge of terrestrial Hg to ocean margins this is what
% reaches the open ocean
k_O_riv_f      = (Te_riv_ocean*fveg + (Te_riv_ocean*fsoil*fCfast)) / Rtf;    % fast
k_O_riv_s      = (Te_riv_ocean*fsoil*fCslow) / Rts;                          % slow
k_O_riv_a      = (Te_riv_ocean*fsoil*fCarmored) / Rta;                       % armored
