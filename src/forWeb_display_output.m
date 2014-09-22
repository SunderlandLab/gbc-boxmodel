%==========================================================================
% OBJECTIVE
%   Print model results to command window. This is really helpful for 
%   testing and debugging.  
%
% REVISION HISTORY
%   08 Sep 2014 - HMA - clean up code and comments for public release
%
% Helen Amos, hamos@hsph.harvard.edu
%==========================================================================

if L_EF;
    
    %---------------------------------------
    % Calculate Anthropogenic Enrichment Values
    %---------------------------------------
    % Enrichment relative to pre-anthropogenic steady state value     
    AtmEF_nat = Matm(end) / Ratm_PI ;
    TfEF_nat  = Mtf(end)  / Rtf_PI  ;
    TsEF_nat  = Mts(end)  / Rts_PI  ;
    TaEF_nat  = Mta(end)  / Rta_PI  ;
    OcsEF_nat = Mocs(end) / Rocs_PI ;
    OciEF_nat = Moci(end) / Roci_PI ;
    OcdEF_nat = Mocd(end) / Rocd_PI ;
    
    if Ldisp;
        disp('-------------------------------------------------------------------')
        disp('ENRICHMENT relative to pre-Anthropogenic steady state (present/pre Anthro era) ')
        disp('-------------------------------------------------------------------')
        disp(['Atmospheric Enrichment                               :   ',num2str((AtmEF_nat))])
        disp(['Fast Terrestrial Enrichment                          :   ',num2str((TfEF_nat))])
        disp(['Slow Terrestrial Enrichment                          :   ',num2str((TsEF_nat))])
        disp(['Armored Terrestrial Enrichment                       :   ',num2str((TaEF_nat))])
        disp(['Surface Ocean Enrichment                             :   ',num2str((OcsEF_nat))])
        disp(['Intermediate Ocean Enrichment                        :   ',num2str((OciEF_nat))])
        disp(['Deep Ocean Enrichment                                :   ',num2str((OcdEF_nat))])
        disp(' ')
    end
    
    %     Enrichment relative to 1840
    indx_1840 = (2008-1840)*(1/dt);
    AtmEF =  Matm(end) / Matm(end-indx_1840) ;
    TfEF  =  Mtf(end)  / Mtf(end-indx_1840)  ;
    TsEF  =  Mts(end)   / Mts(end-indx_1840)  ;
    TaEF  =  Mta(end)  / Mta(end-indx_1840)  ;
    OcsEF =  Mocs(end)  / Mocs(end-indx_1840) ;
    OciEF =  Moci(end)  / Moci(end-indx_1840) ;
    OcdEF =  Mocd(end)  / Mocd(end-indx_1840) ;
    
    if Ldisp;
        disp('-------------------------------------------------------------------')
        disp('ENRICHMENT relative to 1840 (present/1840) ')
        disp('-------------------------------------------------------------------')
        disp(['Atmospheric Enrichment                               :   ',num2str((AtmEF))])
        disp(['Fast Terrestrial Enrichment                          :   ',num2str((TfEF))])
        disp(['Slow Terrestrial Enrichment                          :   ',num2str((TsEF))])
        disp(['Armored Terrestrial Enrichment                       :   ',num2str((TaEF))])
        disp(['Surface Ocean Enrichment                             :   ',num2str((OcsEF))])
        disp(['Intermediate Ocean Enrichment                        :   ',num2str((OciEF))])
        disp(['Deep Ocean Enrichment                                :   ',num2str((OcdEF))])
        disp(' ')
    end
    
end

%---------------------------------------
%Display final reservoir sizes 
%---------------------------------------

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL RESERVOIR SIZES (Mg) ')
    disp('-------------------------------------------------------------------')
    disp(['Atmospheric Reservoir                                :   ',num2str((Matm(end)))])
    disp(['Fast Terrestrial Reservoir                           :   ',num2str(round(Mtf(end)))])
    disp(['Slow Terrestrial Reservoir                           :   ',num2str(round(Mts(end)))])
    disp(['Armored Terrestrial Reservoir                        :   ',num2str(round(Mta(end)))])
    disp(['Surface Ocean Reservoir                              :   ',num2str((Mocs(end)))])
    disp(['Intermediate Ocean Reservoir                         :   ',num2str((Moci(end)))])
    disp(['Deep Ocean Reservoir                                 :   ',num2str(round(Mocd(end)))])
    disp(' ')
end

%---------------------------------------
% CACULATE FINAL YEAR MASS BUDGET FLUXES (Mg/year)
%---------------------------------------

% Atmosphere
DtHgII        =   Matm(end)*k_A_tHgII;       
DtHgIIf       =   fdep_tf*DtHgII;    
DtHgIIs       =   fdep_ts*DtHgII;
DtHgIIa       =   fdep_ta*DtHgII;
DtHg0         =   k_A_tHg0*Matm(end);
DoHgII        =   k_A_oHgII*Matm(end);

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL ATMOSPHERIC DEPOSITION RATES (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp(['Terrestrial Deposition HgII                          :   ',num2str(round(DtHgII))]  )
    disp(['Fast Terrestrial HgII Deposition                     :   ',num2str(round(DtHgIIf))] )
    disp(['Slow Terrestrial HgII Deposition                     :   ',num2str(round(DtHgIIs))] )
    disp(['Armored Terrestrial HgII Deposition                  :   ',num2str(round(DtHgIIa))] )
    disp(['Terrestrial Deposition Hg0                           :   ',num2str(round(DtHg0))]   )
    disp(['Ocean Deposition HgII                                :   ',num2str(round(DoHgII))]  )
    disp(' ')
end

% Fast Terrestrial Reservoir
Tfresp        =   k_Te_rf*Mtf(end);
Tfphoto       =   k_Te_p*Mtf(end);
Tfriv         =   k_T_riv_f*Mtf(end);
Texfs         =   k_T_exfs*Mtf(end);
Texfa         =   k_T_exfa*Mtf(end);

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL FAST TERRESTRIAL RATES (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp(['Fast Terrestrial Respiration                         :   ',num2str(round(Tfresp))]  )
    disp(['Fast Terrestrial Photoreduction and Revolatilization :   ',num2str(round(Tfphoto))] )
    disp(['Fast Terrestrial River Inputs                        :   ',num2str(round(Tfriv))]   )
    disp(['Exchange Fast to Slow Terrestrial                    :   ',num2str(round(Texfs))]   )
    disp(['Exchange Fast to Armored Terrestrial                 :   ',num2str(round(Texfa))]   )
    disp(' ')
end

% Slow Terrestrial Reservoir
Tsresp        =   k_Te_rs*Mts(end);
Tsriv         =   k_T_riv_s*Mts(end);
Texsf         =   k_T_exsf*Mts(end);
Texsa         =   k_T_exsa*Mts(end);

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL SLOW TERRESTRIAL RATES (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp(['Slow Terrestrial Respiration                         :   ',num2str(round(Tsresp))]  )
    disp(['Slow Terrestrial River Inputs                        :   ',num2str(round(Tsriv))]   )
    disp(['Slow Terrestrial Exchange to Fast                    :   ',num2str(round(Texsf))]   )
    disp(['Slow Terrestrial Exchange to Armored                 :   ',num2str(round(Texsa))]   )
    disp(' ')
end

%Armored Soil Reservoir
Taresp        =   k_Te_ra*Mta(end);
Tariv         =   k_T_riv_a*Mta(end);
Texaf         =   k_T_exaf*Mta(end);

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL ARMORED TERRESTRIAL RATES (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp(['Armored Terrestrial Respiration                      :   ',num2str(round(Taresp))] )
    disp(['Armored Terrestrial River Inputs                     :   ',num2str(round(Tariv))]  )
    disp(['Armored Terrestrial Exchange to Fast                 :   ',num2str(round(Texaf))]  )
    disp(' ')
end

% Surface Ocean
OcsUp         =   k_A_oHg0*Matm(end);
OcsEv         =   k_Oc_ev(end)*Mocs(end);
OcsSett       =   k_Oc_sp1*Mocs(end);
OcsDet        =   k_Oc_vsi*Mocs(end);

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL SURFACE OCEAN RATES (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp(['Surface Ocean Uptake of Atmospheric Hg(0)            :   ',num2str(round(OcsUp))]  )
    disp(['Surface Ocean Evasion                                :   ',num2str(round(OcsEv))]  )
    disp(['Particle Settling Surface Ocean                      :   ',num2str(round(OcsSett))])
    disp(['Detrainment to subsurface waters                     :   ',num2str(round(OcsDet))])
    disp(' ')
end

% Intermediate Ocean
OciSett       =   k_Oc_sp2*Moci(end);
OciUpSfc      =   k_Oc_vis*Moci(end);
OciIntDeep    =   k_Oc_vid*Moci(end);

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL INTERMEDIATE OCEAN RATES (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp(['Intermediate Ocean Settling                          :   ',num2str(round(OciSett))]   )
    disp(['Intermediate Ocean Surface Upwelling                 :   ',num2str(round(OciUpSfc))]  )
    disp(['Intermediate Ocean Seawater Sinking                  :   ',num2str(round(OciIntDeep))])
    disp(' ')
end

% Deep Ocean
OcdSett       =   k_Oc_sp3*Mocd(end);
OcdUpInt      =   k_Oc_vdi*Mocd(end);

if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('FINAL DEEP OCEAN RATES (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp(['Deep Ocean Sediment Burial                           :   ',num2str(round(OcdSett))])
    disp(['Deep Ocean Upwelling                                 :   ',num2str(round(OcdUpInt))])
    disp(' ')
end



%%
%--------------------------------------------------------------------------
% Atmosphere: balance emissions and deposition
%--------------------------------------------------------------------------

% no anthropogenic emissions in the pre-1450 run
if ~ L_EF
    Anthro = 0;
end

% calculate total biomass burning, total respiration, and net ocean evasion
if L_EF
    bb_total   = k_Te_BBf_2*Mtf(end) + k_Te_BBs_2*Mts(end) + k_Te_BBa_2*Mta(end); 
else
    bb_total   = k_Te_BBf_1*Mtf(end) + k_Te_BBs_1*Mts(end) + k_Te_BBa_1*Mta(end);
end
resp_total = Tfresp + Tsresp + Taresp;
Ocs_net    = OcsEv - OcsUp;

% total emissions/deposition (Mg/yr)
emis_total = Anthro(end) + E_geo + bb_total + Tfphoto + resp_total + Ocs_net;
dep_total  = DtHgII + DtHg0 + DoHgII;


if Ldisp;
    disp('-------------------------------------------------------------------')
    disp('ATMOSPHERE, Balance emissions and deposition at 2008 (Mg/year) ')
    disp('-------------------------------------------------------------------')
    disp('EMISSIONS ')
    disp( ['  Anthropogenic     : ', num2str( round( Anthro(end) ) ) ] )
    disp( ['  Geogenic          : ', num2str( round( E_geo       ) ) ] )
    disp( ['  Biomass burning   : ', num2str( round( bb_total    ) ) ] )
    disp( ['  Photoreduction    : ', num2str( round( Tfphoto     ) ) ] )
    disp( ['  Respiration       : ', num2str( round( resp_total  ) ) ] )
    disp( ['  Net ocean evasion : ', num2str( round( Ocs_net     ) ) ] )
    disp( ['  TOTAL             : ', num2str( round( emis_total  ) ) ] )
    disp('DEPOSITION')
    disp( ['  Hg(II) to land    : ', num2str( round( DtHgII      ) ) ] )
    disp( ['  Hg(0)  to land    : ', num2str( round( DtHg0       ) ) ] )
    disp( ['  Hg(II) to ocean   : ', num2str( round( DoHgII      ) ) ] )
    disp( ['  TOTAL             : ', num2str( round( dep_total   ) ) ] )
end




