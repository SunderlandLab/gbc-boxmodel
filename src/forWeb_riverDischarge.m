%==========================================================================
% OBJECTIVE
%   Calculate particulate and dissolved Hg export to oceans associated with
%   river discharge. Apply historical scaling factors.
%
% REFERENCES
%    Amos, H. M., et al. (2014), Global biogeochemical implications of 
%      mercury discharges from rivers and sediment burial, Environ. Sci.
%      Technol., 48(16), 9514-9522.
%    Dai, A., and K. E. Trenberth (2002), Estimates of freshwater discharge
%      from continents: Latitudinal and seasonal variations, J. 
%      Hydrometeorol., 3(6), 660-687.
%    Dai, A., et al. (2009), Changes in continental freshwater discharge 
%      from 1948 to 2004, Journal of Climate, 22(10), 2773-2792.
%    Emmerton, C. A., et al. (2013), Mercury export to the Arctic Ocean 
%      from the Mackenzie River, Canada, Environ. Sci. Technol., 47(14), 
%      7644-7654.
%    Hare, A., et al. (2008), Contemporary and preindustrial mass budgets 
%      of mercury in the Hudson Bay marine system: The role of sediment 
%      recycling, Sci. Total Environ., 406(1-2), 190-204.
%    Leitch, D. R., et al. (2007), The delivery of mercury to the Beaufort 
%      Sea of the Arctic Ocean by the Mackenzie River, Sci. Total Environ.,
%      373(1), 178-195.
%    Ludwig, W., and J. L. Probst (1998), River sediment discharge to the 
%      oceans: Present-day controls and global budgets, Am. J. Sci., 298(4), 
%      265-295.
%    Ludwig, W., et al. (2011), Islscp ii global river fluxes of carbon and 
%      sediments to the oceans, edited, Oak Ridge Natinal Laboratory 
%      Distributed Active Archive Center, Oak Ridge, Tennessee, USA.
%    Schuster, P. F., et al. (2011), Mercury export from the Yukon River 
%      basin and potential response to a changing climate, Environ. Sci. 
%      Technol., 45(21), 9262-9267.
%
% REVISION HISTORY
%   12 Jul 2013 - hma - v1.0 created
%   07 Jul 2013 - HMA - use regional Hg concentrations (e.g. N America and
%                       Europe for North Atlantic)
%   19 Aug 2013 - hma - update river concentrations based on more data
%   12 Nov 2013 - hma - revise 1970 to 2008 scaling factors and pre-1970
%                       scaling
%   18 Mar 2014 - HMA - make sure Antarctica is excluded from river inputs
%                       to Southern Ocean since there aren't any rivers there
%   17 Jun 2014 - hma - calculate HgD and Q using 1x1 Q data. The 2x2.5
%                       regridded product has small spatial errors. The
%                       2x2.5 product is OK for plotting with HgP since HgP
%                       is so dominant, but you want to use the 1x1 product
%                       of HgD for the analysis. 
%   17 Sep 2014 - HMA - clean up code and comments for public release
%==========================================================================

%%
%==========================================================================
% Historical scaling factors
%==========================================================================

%----------------------------------------------
% Scaling factors (wrt present-day) for riverine Hg inputs
%----------------------------------------------

% Decades for which you want scaling factors
t_SF = [1960 1970 1980 1990 2000 2008];

% Scaling factors 
% Reference: Table 3 in Amos et al. (2014)

if strcmp(Lscale,'best')
    % best estimate
    SF(1,:) = [7.50 7.50 2.61 1.87 1.33 1]; % North Atlantic, North America
    SF(2,:) = [9.58 9.58 4.80 2.39 1.33 1]; % North Atlantic, Europe
    SF(3,:) = [0.12 0.15 0.25 0.41 0.68 1]; % India
    SF(4,:) = [2.50 4.5  2.15 1.36 1.2  1]; % Mediterranean
    SF(5,:) = [4.25 4.25 1.81 1.44 1.17 1]; % North Pacific, North America
    SF(6,:) = [0.38 0.68 0.62 0.82 0.93 1]; % China    
elseif strcmp(Lscale,'low')
    %low estimate
    SF(1,:) = [3.75 3.75 1.31 1    1    1]; % North Atlantic, North America
    SF(2,:) = [4.52 4.52 2.40 1.98 1    1]; % North Atlantic, Europe
    SF(3,:) = [0.03 0.04 0.06 0.10 0.17 1]; % India
    SF(4,:) = [1.25 2    1.5  1    1    1]; % Mediterranean
    SF(5,:) = [2.13 2.13 1    1    1    1]; % North Pacific, North America
    SF(6,:) = [0.19 0.34 0.31 0.41 0.46 1]; % China
elseif strcmp(Lscale,'high')
    % high estimate
    SF(1,:) = [30   30   3.92 2.81 2    1]; % North Atlantic, North America
    SF(2,:) = [14   14   7.2  2.91 2    1]; % North Atlantic, Europe
    SF(3,:) = [0.50 0.60 1.00 1.64 2.72 1]; % India
    SF(4,:) = [3.75 7    2.8  1.75 1.4  1]; % Mediterranean
    SF(5,:) = [6.38 6.38 2.71 2.15 1.71 1]; % North Pacific, North America
    SF(6,:) = [0.57 1    0.93 1    1    1]; % China
else
    message('Lscale must be best, low, or high!')
end


%----------------------------------------------
% pre-1960 scaled to Hannah's water + soil Hg 
% releases
%----------------------------------------------

% decades for which you want scaling factors
t_SF_pre1960 = 1850:10:1950; 
 
% scaling factors from Horowitz 2013, relative to 1960
SF_pre1960 = [0.114    % 1850
              0.151    % 1860
              0.441    % 1870
              0.575    % 1880
              0.745    % 1890
              0.663    % 1900
              0.522    % 1910
              0.367    % 1920
              0.472    % 1930
              0.823    % 1940
              0.521 ]; % 1950

% Initialize
SF_init = zeros(size(SF,1), numel(SF_pre1960));
SF      = [SF_init SF];

% scale pre-1960 EFs to 1960 EF
for j = 1:6
    SF(j,1:numel(t_SF_pre1960)) = SF(j,numel(t_SF_pre1960)+1)*SF_pre1960.';

end

% time vector associated with scaling factors
t_SF = [t_SF_pre1960 t_SF];

%%
%==========================================================================
% Present-day river inputs
%==========================================================================

%--------------------------------------------------------------------------
% Total suspended sediments (Tg/yr)
%
% Source: ORNL DAAC ISLSCP II project
% http://daac.ornl.gov/cgi-bin/dataset_lister_new.pl?p=29
%
% Notes:
%  The four columns of each file contain:
%  1) Longitude of the center of a cell in decimal degrees.
%  2) Latitude of the center of a cell in decimal degrees. 
%  3) Flux of carbon (in TgC/yr) or sediment (in Tg/yr) into the oceans
%  4) The continental area that is drained to each grid cell
%
% longitude starts at -178.75, increment of 2.5 degrees
% latitude starts  at -77, increment of 2 degrees
%--------------------------------------------------------------------------

% read in data
load ocean_flux_tss_2d.dat

% re-name for ease of coding
COORDraw = ocean_flux_tss_2d(:,1:2);  % [lon,lat]
TSSraw   = ocean_flux_tss_2d(:,3);    % TSS (Tg /yr)

% TSS is currently a 1D vector, I want to make it a 2D matrix

% make vectors of latitude and longitude, grid centers
lat = -89:2:89;
lon = -178.75:2.5:178.75;

% initialize matrix
TSS = zeros( numel(lat), numel(lon) );

% fill in POC matrix
for j = 1:numel(TSSraw)
    
    % for safety's sake
    clear rawlat rawlon;
    
    % lat, lon, and POC value from raw file
    rawlat = COORDraw(j,2); 
    rawlon = COORDraw(j,1); 
    
    % store TSS
    TSS( lat == rawlat, lon == rawlon ) = TSSraw(j);

end

% reshape so that the matrix starts at -180 W, 90 N for plotting
TSS = flipud( TSS );
lat = fliplr( lat );


%%
%--------------------------------------------------------------------------
% Ocean basin mask
%
% Source:
% Created by Helen M Amos
%
% Notes
% (01) Mask was orginally created at 1x1 and regridded to 2x2.5 using the
%      CDO command remapnn.
% (02) Use the 2x2.5 mask for TSS/HgP and the 1x1 mask for Q/HgD. (17 Jun
%      2014, hma) 
%
% Ocean basins:
%    land = 1
%    ARC : Arctic Ocean     = 2
%    NPA : North Pacific    = 3
%    NAT : North Atlantic   = 4
%    SOC : Southern Ocean   = 5
%    PAC : Pacific Ocean    = 6
%    SAT : South Atlantic   = 7
%    IND : Indian Ocean     = 8
%    MED : Mediterranean    = 9
%--------------------------------------------------------------------------

% open file
clear filename filename1deg;
filename     = 'oceanmask_2x2.5.nc';                    % 2x2.5
filename1deg = 'oceanmask_1x1.nc';                      % 1x1
ncid         = netcdf.open( filename, 'NOWRITE' );      % 2x2.5
ncid1deg     = netcdf.open( filename1deg, 'NOWRITE' );  % 1x1

% read in ocean mask
mask         = ncread( filename     , 'ocean mask' );   % 2x2.5
mask1deg     = ncread( filename1deg , 'ocean mask' );   % 1x1

% close file
netcdf.close( ncid );
netcdf.close( ncid1deg );

% for safety's sake
clear ncid ncid1deg;

% set land to 0
mask( mask < 2 )         = 0; % 2x2.5
mask1deg( mask1deg < 2 ) = 0; % 1x1

% re-shape so that matrix starts at -180 W, 90N
mask     = mask.';
mask     = flipud( mask );
mask1deg = mask1deg.';
mask1deg = flipud( mask1deg );


% Push the ocean basin mask inland so you don't accidentally exclude any
% coastal boxes.

% force land box to be an ocean box
nrow = size(mask,1);  % number of rows
ncol = size(mask,2);  % number of columns

% dummy matrix
dummy4 = ones( size(mask,1), size(mask,2)) ;

% look at all the cells and record the location of coastal cells and save
% the ID number of the adjacent basin

% loop over rows
for k = 1:ncol;   
for j = 2:nrow-1; 
    if mask(j,k) == 0 && mask(j-1,k) ~= 0 && mask(j+1,k) ==0
        dummy4(j,k) = mask(j-1,k);
    end

    if mask(j,k) == 0 && mask(j-1,k) == 0 && mask(j+1,k) ~=0
        dummy4(j,k) = mask(j+1,k);
    end
end
end

% loop over columns
for k = 2:ncol-1;
for j = 1:nrow;
    if mask(j,k) == 0 && mask(j,k-1) ~= 0 && mask(j,k+1) ==0
        dummy4(j,k) = mask(j,k-1);
    end

    if mask(j,k) == 0 && mask(j,k-1) == 0 && mask(j,k+1) ~=0
        dummy4(j,k) = mask(j,k+1);
    end

end
end


% force land to be 1
mask( mask == 0 ) = 1;

% include coastline in basins
mask = mask.* dummy4;


% display global and basin  TSS totals
if Ldisp;
    disp(' ')
    disp(' ')
    disp('TSS totals (Tg /yr)')
    disp(['   global = ', num2str( sum(sum( TSS )) )])
    basinNames = {'ARC','NPA','NAT','SOC','PAC','SAT','IND','MED'};
    for j = 2:9
        disp(['   ',basinNames{j-1},'    = ', num2str( sum(sum( TSS(mask == j) )) )])
    end
    disp(' ')
    disp(' ')
end


%%
%--------------------------------------------------------------------------
% Annual mean river discharge (Sv)
%
% Data source:
% Dai & Trenberth (2002), Dai et al. (2009)
%
% URL:
% http://www.cgd.ucar.edu/cas/catalog/surface/dai-runoff/index.html
%
% Download date: 05 April 2013
%
% Notes:
%   (01) Original file is in unformatted FORTRAN binary. Lee T. Murray
%        wrote an IDL script to read in the binary file and I wrote an IDL
%        script to save it to netCDF. 
%
%--------------------------------------------------------------------------

% Add in-land freshwater discharge to coastal boxes. 
%
% Note from Dai & Trenberth 2009, on their data website:
% "The values are discharge in Sv or 1E6 m3/s) within the coastal boxes 
% (some boxes may be at the farthest downstream stations of small rivers 
% and thus not right on the coasts in the 921 river case. These boxes need 
% to be included in computing discharge from a lat zone)."
%
% Do the adding by saving out the data to a .txt file, then opening it in 
% Excel and modifying the file manually. Add inland boxes to the nearest
% coastal box w/ a non-zero Q at the same latitude. Or if the closest
% non-zero Q is an adjacent box in the same longitude, add it there. 
%
load 'river_Q_1x1_ADDinlandQ_20140617.txt'                        % after adding inland boxes
river_ann1deg = river_Q_1x1_ADDinlandQ_20140617;

% Extend 1x1 mask in land so that it overlaps with all Q points

% force land box to be an ocean box
nrow1deg = size(mask1deg,1);  % number of rows
ncol1deg = size(mask1deg,2);  % number of columns

% look at all the cells and record the location of coastal cells and save
% the ID number of the adjacent basin
%
% inch the ocean mask inland over 3 iterations
for s=1:3;

    % dummy matrix
    clear dummy4; % for safety's sake
    dummy4 = ones( size(mask1deg,1), size(mask1deg,2)) ;

    
    for k = 1:ncol1deg;
        for j = 2:nrow1deg-1;
            if mask1deg(j,k) == 0 && mask1deg(j-1,k) ~= 0 && mask1deg(j+1,k) ==0
                dummy4(j,k) = mask1deg(j-1,k);
            end
            
            if mask1deg(j,k) == 0 && mask1deg(j-1,k) == 0 && mask1deg(j+1,k) ~=0
                dummy4(j,k) = mask1deg(j+1,k);
            end
        end
    end
    
    %loop over columns
    for k = 2:ncol1deg-1;
        for j = 1:nrow1deg;
            if mask1deg(j,k) == 0 && mask1deg(j,k-1) ~= 0 && mask1deg(j,k+1) ==0
                dummy4(j,k) = mask1deg(j,k-1);
            end
            
            if mask1deg(j,k) == 0 && mask1deg(j,k-1) == 0 && mask1deg(j,k+1) ~=0
                dummy4(j,k) = mask1deg(j,k+1);
            end
            
        end
    end
    
    % force land to be 1
    mask1deg( mask1deg == 0 ) = 1;
    
    % include coastline in basins
    mask1deg = mask1deg.* dummy4;
    
    % force land back to 0
    mask1deg( mask1deg == 1 ) = 0;
    
    
end;

% force land to be 1
mask1deg( mask1deg == 0 ) = 1;


%%
%==========================================================================
% Estimate HgD and HgP (Mg/yr) based on observations 
%==========================================================================

% Make a matrix of repeating longitudes and latitudes
% 
% 2x2.5 resolution
lon_mat = lon;
for j = 1:size(TSS,1)-1
    lon_mat = vertcat( lon_mat, lon);
end

lat_mat = lat.';
for j = 1:size(TSS,2)-1
    lat_mat = horzcat( lat_mat, lat.');
end
%
% 1x1 resolution
lon1deg     = [-179.5:179.5];     % vector of longitudes at grid box centers
lon_mat1deg = lon1deg;            % intialize
for j = 1:size(mask1deg,1)-1
    lon_mat1deg = vertcat( lon_mat1deg, lon1deg);
end

lat1deg     = [89.5:-1:-89.5].';  % vector of latitudes at grid box centers
lat_mat1deg = lat1deg;            % intialize    
for j = 1:size(mask1deg,2)-1
    lat_mat1deg = horzcat( lat_mat1deg, lat1deg);
end

% Use Q (km3/yr) to convert HgD (pM) to Mg/yr:
%
%  pmol     km3     1e-12 mol     1e12 L     201 g     1 Mg
% ------ x ----- x ----------- x -------- x ------- x ------
%   L        yr      1 pmol       1 km3      1 mol     1e6 g
%
% ==> HgD * Q * 201 *1e-6
%
pM_MgYr = 201 * 1e-6;

% Use TSS (g /yr) to convert HgP (pmol/g) to Mg/yr:
%
%  nmol     g       1e-9 mol     201 g     1 Mg     
% ------ x ----- x ---------- x ------- x ------ 
%   g       yr        nmol        mol      1e6 g
%
% ==> HgP * TSS * 201 * 1e-15
%
nmolg_MgYr = 201 * 1e-3;


%--------------------------------------------------------------------------
% Flow-weighted mean estimates of HgD (pM) and HgP (nmol/g) based on 
% published literature. These represent PRESENT-DAY concentratations
%--------------------------------------------------------------------------

% recorded here as [mean, lower uncertainty bound, upper uncertainty bound],
% with the full range of the observations noted in the comments

if strcmp(Lriver,'best') 
    riv_idx = 1;
elseif strcmp(Lriver,'low')
    riv_idx = 2;
elseif strcmp(Lriver,'high')
    riv_idx = 3;
else
    message('Options are: high, low, best')
end

%---------------------------
% Mediterranean
%---------------------------
LHgP_nmolg_MED = [0.10, 0.10-0.09, 0.10+0.09];            
LHgD_pM_MED    = [1.8 , 1.8-0.8  , 1.8+0.8  ];          

%---------------------------
% Indian Ocean: India, Bangladesh, and Myanmar
%---------------------------
LHgP_nmolg_IND1 = [2.69, 2.69-2.69, 2.69+2.86];    
LHgD_pM_IND1    = [50  , 50-33    , 50+33    ]; 

%---------------------------
% Arctic Ocean
%---------------------------

% North America and Europe 
LHgP_nmolg_ARC1 = [0.37, 0.37-0.32, 0.37+0.32];      
LHgD_pM_ARC1    = [6.8 , 6.8-1.7  , 6.8+1.7  ];

% Russia
LHgP_nmolg_ARC2 = [0.46, 0.46-0.38, 0.46+0.38]; 
LHgD_pM_ARC2    = [8.5 , 8.5-0.9  , 8.5+0.9  ];

%---------------------------
% North Atlantic
%---------------------------

% North America
LHgP_nmolg_NAT1 = [0.47, 0.47-0.32, 0.47+0.32];
LHgD_pM_NAT1    = [8.7 , 8.7-5.7  , 8.7+5.7  ];

% Europe
LHgP_nmolg_NAT2 = [0.49, 0.49-0.42, 0.49 + 0.42];
LHgD_pM_NAT2    = [9.1 , 9.1-2.1  , 9.1+2.1    ];

%---------------------------
% South Atlantic
%---------------------------

% North America (based on the Mississippi)
LHgP_nmolg_SAT1 = [0.30 ,0.30-0.24, 0.30+0.24 ];       
LHgD_pM_SAT1    = [5.5  , 5.5-3.6, 5.5+3.6];

% South America 
LHgP_nmolg_SAT2 = [1.50, 1.50-1.24, 1.50+1.24];
LHgD_pM_SAT2    = [28  , 28-18    , 28+18    ];

% Africa (based on the Nile)
LHgP_nmolg_SAT3 = [0.11, 0.11-0.10, 0.11+0.10];
LHgD_pM_SAT3    = [2   , 2-1.3    , 2+1.3    ];


%---------------------------
% North Pacific
%---------------------------

% North America
LHgP_nmolg_NPA1 = [0.34, 0.34-0.33, 0.34+0.33];         
LHgD_pM_NPA1    = [6.4 , 6.4-3.2  , 6.4+3.2  ];         

% Asia
LHgP_nmolg_NPA2 = [5.91, 5.91-5.76, 5.91+5.76];
LHgD_pM_NPA2    = [110 , 110-55   , 110+55   ];              


%---------------------------
% Pacifc
%---------------------------

% SE Asia (based on obs from Mekong River, Vietnam)
LHgP_nmolg_PAC1 = [0.19, 0.19-0.19, 0.19+0.19];
LHgD_pM_PAC1    = [3.6 , 3.6-2.3  , 3.6+2.3  ];    


%--------------------------------------------------------------------------
% Calculate Hg inputs from rivers to each ocean basin
%
% Scale present-day osbervations to estimate historical Hg concentrations 
%--------------------------------------------------------------------------

% initialize
river_HgD_MgYr_save  = zeros( 8, numel(t_SF) );     % HgD inputs to each basin (Mg/yr)
river_HgP_MgYr_save  = zeros( 8, numel(t_SF) );     % HgD inputs to each basin (Mg/yr)


for j = 1:numel(t_SF); % loop over decades for scaling
    
    % for safety's sake
    river_HgD       = 0 * river_ann1deg; % in Mg/yr @ 1x1
    river_HgP       = 0 * TSS;     % in Mg/yr @ 2x2.5
    
    %----------------------
    % Mediterranean
    river_HgP( mask == 9 )     = SF(4,j) * LHgP_nmolg_MED(riv_idx) * nmolg_MgYr * TSS( mask == 9 );
    river_HgD( mask1deg == 9 ) = SF(4,j) * LHgD_pM_MED(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 9 );
       
    %----------------------
    % Indian Ocean
    
    % India, Bangladesh, and Myanmar
    river_HgP( mask == 8 & lon_mat > 47 & lat_mat > 5 )     = SF(3,j) * LHgP_nmolg_IND1(riv_idx) * nmolg_MgYr * TSS( mask == 8 & lon_mat > 47 & lat_mat > 5 );
    river_HgD( mask1deg == 8 & lon_mat1deg > 47 & lat_mat1deg > 5 ) = SF(3,j) * LHgD_pM_IND1(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 8 & lon_mat1deg > 47 & lat_mat1deg > 5 );
    
    % SE Asia
    %   No historical information available, so no scaling applied
    river_HgP( mask == 8 & lon_mat > 80 & lat_mat > -10 )     = LHgP_nmolg_PAC1(riv_idx) * nmolg_MgYr * TSS( mask == 8 & lon_mat > 80 & lat_mat > -10 );
    river_HgD( mask1deg == 8 & lon_mat1deg > 80 & lat_mat1deg > -10 ) = LHgD_pM_PAC1(riv_idx)  * pM_MgYr * river_ann1deg( mask1deg == 8 & lon_mat1deg > 80 & lat_mat1deg > -10 );
    
    % Australia
    %   No obs, so assume same value as SE Asia
    %   No historical information available, so no scaling applied
    river_HgP( mask == 8 & lon_mat > 80 & lat_mat <= -10 )     = LHgP_nmolg_PAC1(riv_idx) * nmolg_MgYr * TSS( mask == 8 & lon_mat > 80 & lat_mat <= -10 );
    river_HgD( mask1deg == 8 & lon_mat1deg > 80 & lat_mat1deg <= -10 ) = LHgD_pM_PAC1(riv_idx)  * pM_MgYr * river_ann1deg( mask1deg == 8 & lon_mat1deg > 80 & lat_mat1deg <= -10 );
    
    % Africa
    %   No historical information available, so no scaling applied
    river_HgP( mask == 8 & lon_mat <= 80 & lat_mat <= 5 )     = LHgP_nmolg_SAT3(riv_idx) * nmolg_MgYr * TSS( mask == 8 & lon_mat <= 80 & lat_mat <= 5 );
    river_HgD( mask1deg == 8 & lon_mat1deg <= 80 & lat_mat1deg <= 5 ) = LHgD_pM_SAT3(riv_idx)   * pM_MgYr * river_ann1deg( mask1deg == 8 & lon_mat1deg <= 80 & lat_mat1deg <= 5 );
   
    
    %----------------------                     
    % Arctic
    
    % North America and Europe (lon < 30 E)
    % No historical information available, so no scaling applied
    river_HgP( mask == 2 & lon_mat < 30 )     = LHgP_nmolg_ARC1(riv_idx) * nmolg_MgYr * TSS( mask == 2 & lon_mat < 30 );
    river_HgD( mask1deg == 2 & lon_mat1deg < 30 ) = LHgD_pM_ARC1(riv_idx)   * pM_MgYr * river_ann1deg( mask1deg == 2 & lon_mat1deg < 30);
    
    % Russia
    %   No historical information on Arctic
    river_HgP( mask == 2 & lon_mat >= 30 )     = LHgP_nmolg_ARC2(riv_idx) * nmolg_MgYr * TSS( mask == 2 & lon_mat >= 30 );
    river_HgD( mask1deg == 2 & lon_mat1deg >= 30 ) = LHgD_pM_ARC2(riv_idx)   * pM_MgYr * river_ann1deg( mask1deg == 2 & lon_mat1deg >= 30);
    
   
    %----------------------
    % North Atlantic
    
    % North America
    river_HgP( mask == 4 & lon_mat < -52.5)     = SF(1,j) * LHgP_nmolg_NAT1(riv_idx) * nmolg_MgYr * TSS( mask == 4 & lon_mat < -52.5 );
    river_HgD( mask1deg == 4 & lon_mat1deg < -52.5) = SF(1,j) * LHgD_pM_NAT1(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 4 & lon_mat1deg < -52.5 );
    
    % Europe
    river_HgP( mask == 4 & lon_mat >= -52.5)     = SF(2,j) * LHgP_nmolg_NAT2(riv_idx) * nmolg_MgYr * TSS( mask == 4 & lon_mat >= -52.5 );
    river_HgD( mask1deg == 4 & lon_mat1deg >= -52.5) = SF(2,j) * LHgD_pM_NAT2(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 4 & lon_mat1deg >= -52.5 );
 

    %---------------------- 
    % South Atlantic

    % North America
    %   Same scaling as North America to North Atlantic
    river_HgP( mask == 7 & lon_mat < -78 & lat_mat > 24 )     = SF(1,j) * LHgP_nmolg_SAT1(riv_idx) * nmolg_MgYr * TSS( mask == 7 & lon_mat < -78 & lat_mat > 24 );
    river_HgD( mask1deg == 7 & lon_mat1deg < -78 & lat_mat1deg > 24 ) = SF(1,j) * LHgD_pM_SAT1(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 7 & lon_mat1deg < -78 & lat_mat1deg > 24 );

    % South America and Africa
    % No historical information available, so no scaling applied
    %
    %-- NW Africa
    river_HgP( mask == 7 & lon_mat >= -78 & lat_mat > 24 ) = LHgP_nmolg_SAT3(riv_idx) * nmolg_MgYr * TSS( mask == 7 & lon_mat >= -78 & lat_mat > 24 );
    river_HgD( mask1deg == 7 & lon_mat1deg >= -78 & lat_mat1deg > 24 ) = LHgD_pM_SAT3(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 7 & lon_mat1deg >= -78 & lat_mat1deg > 24 );

    %-- S America 
    river_HgP( mask == 7 & lon_mat < -20 & lat_mat <= 24 )     = LHgP_nmolg_SAT2(riv_idx) * nmolg_MgYr * TSS( mask == 7 & lon_mat < -20 & lat_mat <= 24 );
    river_HgD( mask1deg == 7 & lon_mat1deg < -20 & lat_mat1deg <= 24 ) = LHgD_pM_SAT2(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 7 & lon_mat1deg < -20 & lat_mat1deg <= 24 );

    %-- the rest of Africa
    river_HgP( mask == 7 & lon_mat >= -20 & lat_mat <= 24 )     = LHgP_nmolg_SAT3(riv_idx) * nmolg_MgYr * TSS( mask == 7 & lon_mat >= -20 & lat_mat <= 24 );
    river_HgD( mask1deg == 7 & lon_mat1deg >= -20 & lat_mat1deg <= 24 ) = LHgD_pM_SAT3(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 7 & lon_mat1deg >= -20 & lat_mat1deg <= 24 );
    
                         
    %---------------------- 
    % North Pacific
    
    % North America
    river_HgP( mask == 3 & lon_mat <= -90 )     = SF(5,j) * LHgP_nmolg_NPA1(riv_idx) * nmolg_MgYr * TSS( mask == 3 & lon_mat <= -90 );
    river_HgD( mask1deg == 3 & lon_mat1deg <= -90 ) = SF(5,j) * LHgD_pM_NPA1(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 3 & lon_mat1deg <= -90 );

    % Asia
    %   Use the scale factors for China.
    river_HgP( mask == 3 & lon_mat >= 90 )     = SF(6,j) * LHgP_nmolg_NPA2(riv_idx)  * nmolg_MgYr * TSS( mask == 3 & lon_mat >= 90 );
    river_HgD( mask1deg == 3 & lon_mat1deg >= 90 ) = SF(6,j) * LHgD_pM_NPA2(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 3 & lon_mat1deg >= 90 );

    
    %---------------------- 
    % Pacifc
    
    % China
    %   Same HgD and HgP concentrations as Chinese rivers entering the
    %   North Pacific. 
    river_HgP( mask == 6 & lon_mat > 0 & lat_mat >= 15 )     = SF(6,j) *  LHgP_nmolg_NPA2(riv_idx) * nmolg_MgYr * TSS( mask == 6 & lon_mat > 0 & lat_mat >= 15);
    river_HgD( mask1deg == 6 & lon_mat1deg > 0 & lat_mat1deg >= 15 ) = SF(6,j) *  LHgD_pM_NPA2(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 6 & lon_mat1deg > 0 & lat_mat1deg >= 15);
  
    % SE Asia (based on obs from Mekong River, Vietnam)
    %   No historical information available, so no scaling applied
    river_HgP( mask == 6 & lon_mat > 0 & lat_mat > -10 & lat_mat < 15 )     =  LHgP_nmolg_PAC1(riv_idx) * nmolg_MgYr * TSS( mask == 6 & lon_mat > 0 & lat_mat > -10 & lat_mat < 15);
    river_HgD( mask1deg == 6 & lon_mat1deg > 0 & lat_mat1deg > -10 & lat_mat1deg < 15 ) =  LHgD_pM_PAC1(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 6 & lon_mat1deg > 0 & lat_mat1deg > -10 & lat_mat1deg < 15);

    
    % Australia
    %   No HgD or HgP obs, assume to be the same as SE Asia
    river_HgP( mask == 6 & lon_mat > 0 & lat_mat <= -10 )     =  LHgP_nmolg_PAC1(riv_idx) * nmolg_MgYr * TSS( mask == 6 & lon_mat > 0 & lat_mat <= -10 );
    river_HgD( mask1deg == 6 & lon_mat1deg > 0 & lat_mat1deg <= -10 ) =  LHgD_pM_PAC1(riv_idx) * pM_MgYr* river_ann1deg( mask1deg == 6 & lon_mat1deg > 0 & lat_mat1deg <= -10 );
 
    % South America
    %   Assume to be same as S America to S. Atlantic
    river_HgP( mask == 6 & lon_mat < 0 )     = LHgP_nmolg_SAT2(riv_idx) * nmolg_MgYr * TSS( mask == 6 & lon_mat < 0 );
    river_HgD( mask1deg == 6 & lon_mat1deg < 0 ) = LHgD_pM_SAT2(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 6 & lon_mat1deg < 0 );

    
    %---------------------- 
    % Southern Ocean
    %   - Assume the same HgD and HgP values as S. America, since that's most of the TSS and Q source.
    %   - make sure Antarctica is exclude, since there aren't any rivers there
    %   - No historical information available, so no scaling applied
    river_HgP( mask == 5 & lat_mat > -60) = LHgP_nmolg_SAT2(riv_idx) * nmolg_MgYr * TSS( mask == 5 & lat_mat > -60);
    river_HgD( mask1deg == 5 & lat_mat1deg > -60) = LHgD_pM_SAT2(riv_idx) * pM_MgYr * river_ann1deg( mask1deg == 5 & lat_mat1deg > -60);
    

    %----------------------------------------------------
    % display global and basin river HgD and HgP totals
    %----------------------------------------------------
    if Ldisp;
        disp(' ')
        disp(['Estimate: ',Lriver] )
        disp(['DECADE: ', num2str(t_SF(j)) ])
        disp(' ')
        disp(' ')
        disp('           HgD in Mmol/yr (in Mg/yr)  | HgP  in Mmol/yr (in Mg/yr)  ')
        
        for m = 2:9
            disp(['   ',basinNames{m-1},'   '...
                num2str( sum(sum( river_HgD(mask1deg == m) ))/200.59 ),'     (',...
                num2str( sum(sum( river_HgD(mask1deg == m) )) ),')   |  '      ,...
                num2str( sum(sum( river_HgP(mask == m) ))/200.59 ) ,'   ('     ,...
                num2str( sum(sum( river_HgP(mask == m) )) ),')' ])
        end
        disp(['   global : HgD = ', num2str( sum(sum( river_HgD )) ), ' (Mg/yr),   HgP = ',...
            num2str( sum(sum( river_HgP )) ), ' (Mg/yr)  ' ])
        
        disp(' ')
        disp(' ')
    end
    
    % loop over ocean basins
    for w = 2:9
        % HgD input to each basin basin (Mg/yr)
        river_HgD_MgYr_save(w-1,j) = sum( sum( river_HgD(mask == w) ) ); 

        % HgP input to each basin basin (Mg/yr)
        river_HgP_MgYr_save(w-1,j) = sum( sum( river_HgP(mask == w) ) ); 
    end
end


%%
%--------------------------------------------------------------------------
% Estimate "background"  river Hg inputs (i.e. what's coming out of rivers
% in the absence of direct dumping of industrial and municipal wastewaters,
% so just reflecting natural Hg and enrichment from atmospheric
% deposition).
%--------------------------------------------------------------------------

% Yukon River Hg 
% Ref: Schuster et al., 2011, ES&T
HgD_yukon = 1.9;    % ng/L
HgP_yukon = 13.1;   % ng/L
Q_yukon   = 203;    % km3/yr
TSS_yukon = 6.1e10; % kg/yr

% convert ng/L to pM and nmol/g
HgD_yukon = HgD_yukon*1e3/201;  % pM
HgP_yukon = HgP_yukon*1e9*Q_yukon/TSS_yukon/201;

% Mackenzie River Hg
% Ref: Leitch et al., 2007
HgD_macken1 = 14;     % 14 +/- 11 pM
HgP_macken1 = 0.05;   % 0.05 + 0.5 nmol/g
Q_macken1   = 330;    % km3/yr
TSS_macken1 = 1.3e11; % kg/yr  ref: Hare et al. 2008

% Mackenzie River Hg
% ref: Emmerton et al., 2013, ES&T
HgD_macken2 = 8.0;              % 8.0 +/- 3.0 pM
HgP_macken2 = (0.6+0.1)/2;      % mean of range (0.1 - 0.6) nmol/g. reported as 88 +/- 5% of THg (in ng/L)
Q_macken2   = 290-330;          % km3/yr
TSS_macken2 = (160e9 + 40e9)/2; % mean of range (40 - 160e9) kg/yr

% weighted mean of Yukon and Mackenzie [Hg]
CHgD_pristine = ( HgD_yukon*Q_yukon + HgD_macken1*Q_macken1 + HgD_macken2*Q_macken1 )...
               /(Q_yukon + 2*Q_macken1); % pM
           
CHgP_pristine = ( HgP_yukon*Q_yukon + HgP_macken1*Q_macken1 + HgP_macken2*Q_macken1 )...
               /(Q_yukon + 2*Q_macken1); % nmol/g
           

% estimate global background Hg inputs from rivers (Mg/yr)
% -- the assumption here is that the Hg concentrations in pristine Arctic
%    rivers are representative of background [Hg] everywhere
IHgD_pristine = CHgD_pristine * sum(sum(river_ann1deg)) * pM_MgYr;
IHgP_pristine = CHgP_pristine * sum(sum(TSS      )) * nmolg_MgYr;

ITHg_pristine = IHgD_pristine + IHgP_pristine;

