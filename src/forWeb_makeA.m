%==========================================================================
% OBJECTIVE
%   Assemble matrix A, such that dM/dt = A*M + E
%
% REVISION HISTORY
%   18 Jan 2012 - hma - version 1.0 written
%   19 Jul 2012 - hma - remove deep mineral pool and treat geogenic
%                       emissions as an external forcing
%   20 Mar 2014 - hma - update parameterization for rivers and add burial
%                       in coastal benthic sediment
%   08 Sep 2014 - HMA - clean up code and comments for public release
%   08 Sep 2014 - HMA - clean up code and comments for public release
%
% Helen M. Amos, hamos@hsph.harvard.edu
%==========================================================================

% clear variable, for safety's sake
clear A;

% First-order rate coefficients, k (1/yr)
forWeb_rate_coeffs

% The rate coefficient for biomass burning differs between the
% anthropogenic and pre-anthropogenic simulations. It's the only rate
% coefficient with a time-dependence. 
if sim_type == 1;     % <-- pre-anthropogenic era
    k_Te_BBf = k_Te_BBf_1; % fast soil
    k_Te_BBs = k_Te_BBs_1; % slow soil
    k_Te_BBa = k_Te_BBa_1; % armored soil
    
elseif sim_type == 2; % <-- anthropogenic era
    k_Te_BBf = k_Te_BBf_2; % fast soil
    k_Te_BBs = k_Te_BBs_2; % slow soil
    k_Te_BBa = k_Te_BBa_2; % armored soil    
else
    message('Invalid simulation type! Must be 1 or 2.')
end

%--------------------------------------------------------------------------
% Assemble matrix A, such that dM/dt = A*M + E
%--------------------------------------------------------------------------

%-- atmosphere
Aatm = [-(k_A_oHgII + k_A_oHg0 + k_A_tHgII + k_A_tHg0)   % Matm term
         (k_Te_rf + k_Te_p + k_Te_BBf)                   % Mtf term
         (k_Te_rs + k_Te_BBs)                            % Mts
         (k_Te_ra + k_Te_BBa)                            % Mta
          k_Oc_ev                                        % Mocs
          0                                              % Moci
          0 ];                                           % Mocd
      
%-- fast soil pool
Atf = [ (k_A_tHgII * fdep_tf + k_A_tHg0)
       -(k_T_riv_f + k_Te_rf + k_Te_p + k_T_exfs + k_T_exfa + k_Te_BBf);
         k_T_exsf
         k_T_exaf
         0
         0
         0 ];
     
%-- slow soil pool
Ats = [ k_A_tHgII*fdep_ts
        k_T_exfs 
      -(k_Te_rs + k_T_exsf + k_T_exsa + k_T_riv_s + k_Te_BBs)
        0
        0
        0
        0 ];
    
%-- amored soil pool
Ata = [ k_A_tHgII * fdep_ta
        k_T_exfa
        k_T_exsa
      -(k_Te_ra + k_T_exaf + k_T_riv_a + k_Te_BBa)
        0
        0
        0 ];
    
%-- surface ocean
Aocs = [ k_A_oHgII + k_A_oHg0
         k_O_riv_f
         k_O_riv_s
         k_O_riv_a
       -(k_Oc_sp1 + k_Oc_ev + k_Oc_vsi)
         k_Oc_vis
         0];
     
%-- subsurface ocean
Aoci = [ 0
         0
         0
         0
         k_Oc_sp1 + k_Oc_vsi
       -(k_Oc_vis+ k_Oc_vid + k_Oc_sp2)
         k_Oc_vdi ];
     
%-- deep ocean
Aocd = [ 0
         0
         0
         0
         0
         k_Oc_vid + k_Oc_sp2
       -(k_Oc_vdi + k_Oc_sp3)];
     
     
%-- matrix A     
A = [Aatm.' ; Atf.' ; Ats.' ; Ata.' ; Aocs.' ; Aoci.' ; Aocd.'];
