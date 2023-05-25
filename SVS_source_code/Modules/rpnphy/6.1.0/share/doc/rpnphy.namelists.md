### convection_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| bkf_closures | Select closures for shallow convection<br>- 'CAPE'<br>- 'EQUILIBRIUM' | 'CAPE' | character(len=16) |
| bkf_detrains | Select formulation of fractional detrainment rate for shallow convection<br>- 'BECHTOLD01'<br>- 'CUIJPERS95'<br>- 'DEROOY10' | 'BECHTOLD01' | character(len=16) |
| bkf_entrains | Select formulation of fractional entrainment rate for shallow convection<br>- 'BECHTOLD01'<br>- 'BECHTOLD08'<br>- 'DEROOY11'<br>- 'SIEBESMA03' | 'BECHTOLD01' | character(len=16) |
| bkf_evaps | Evaporate detrained condensate in shallow convection | .false. | logical |
| bkf_kch | Number of species for convective transport (never tested) | 0 | integer |
| bkf_kens | Number of additional ensemble members (max 3) for deep bkf convection | 0 | integer |
| bkf_kice | Take ice phase into account in deep bkf (yes=1) | 1 | integer |
| bkf_ktdia | Limit vertical computation by ktdia-1 levels | 1 | integer |
| bkf_lch1conv | Activate convective transport of species for deep and shallow bkf | .false. | logical |
| bkf_ldown | Allow downdrafts in deep bkf | .true. | logical |
| bkf_lshalm | Activate shallow convective momentum transport | .false. | logical |
| bkf_rads | Cloud radii at LCL for bkf_shallow<br>from bkf_rads(1) to bkf_rads(2)<br>with increment bkf_rads(3) | (/50., 50., 0./) | real |
| bkf_tperts | Temperature perturbations at LCL for triggering bkf_shallow<br>An ensemble of shall. cumuli will be generated using perturbations<br>starting from  bkf_tperts(1) to bkf_tperts(2)<br>with increment bkf_tperts(3) | (/0.2, 0.2, 0./) | real |
| deep | Deep convection scheme name<br>- 'NIL     ' :<br>- 'SEC     ' :<br>- 'OLDKUO  ' :<br>- 'KUOSTD  ' :<br>- 'KFC     ' :<br>- 'KFC2    ' :<br>- 'KFC3    ' :<br>- 'BECHTOLD' : | 'nil' | character(len=16) |
| deep_cloudobj | Treat convective clouds as cloud objects | .false. | logical |
| deep_codecay | Decay timescale for convective cloud objects (seconds) | 600. | real |
| deep_conserve | Conservation corrections for deep convective scheme<br>- 'NIL   ' : No conservation correction applied<br>- 'TEND  ' : Temperature and moisture tendencies corrected<br>- 'PRECIP' : Surface precipitation rate corrected | 'PRECIP' | character(len=16) |
| deep_timeconv |  | -1. | real |
| deep_timeent |  | -1. | real |
| deep_timerefresh |  | -1. | real |
| kfcdepth | Minimum depth of conv. updraft for KFC  trigger (m) | 4000. | real |
| kfcdpdd | Maximum depth of the downdraft detrainment layer (Pa) for 'kfc2' and 'kfc3' | 10000. | real |
| kfcmom | generate wind tendencies in KFC or deep BKF if .true. | .false. | logical |
| kfcprod | Compute production terms for Kain-Fritsch scheme | .false. | logical |
| kfcrad | Initial convective updraft radius in KFC scheme(m) | 1500. | real |
| kfcradw | Convective updraft radius over water in KFC scheme(m) | -1. | real |
| kfctaucape | Varies convective timescale as a function of CAPE for Kain-Fritsch scheme<br>KFCTAUCAPE = time1, time2, cmean, dcape<br>- time1 (s): max kfctimec<br>- time2 (s): min kfctimec<br>- cmean (J/Kg): cape value at which kfctimec will be mean of time1 and time2<br>- dcape (J/Kg): decrease in kfctimec from time1 to time2 will occur over range cmean-dcape to cmean+dcape | (/-1., -1., -1., -1./) | real |
| kfctrig4 | Trigger parameter of Kain-Fritsch convection scheme (WKLCL).<br>Trigger parameter will increase from kfctrig4(3) to kfctrig4(4) [m/s]<br>between timestep kfctrig4(1) and timestep kfctrig4(2) | (/0., 0., 0.05, 0.05/) | real |
| kfctriga | Nominal resolution for which KFCTRIG4 is set.<br>This is inactive if value <= 0. | -1.0 | real |
| kfctrigl | Over land and lakes we keep the value set by the "ramp" above over sea water:<br>- for :lat: >= TRIGLAT(2) we keep value set by the "ramp" KFCTRIG4<br>- for :lat: <= TRIGLAT(1) we use the new value KFCTRIGL [m/s]<br>- and linear interpolation in between TRIGLAT(1) and TRIGLAT(2) | 0.05 | real |
| kfctriglat | Logical key for variation of the trigger function depending on latitude and land-sea-lake mask | .false. | logical |
| kfctrigtau | Relaxation timescale for trigger velocity | -1. | real |
| kfctrigw | Trigger parameter of Kain-Fritsch convection scheme (WKLCL).<br>Trigger parameter will increase from kfctrigw(3) to kfctrigw(4) [m/s]<br>between wstar values kfctrigw(1) and kfctrigw(2) | (/0., 0., 0., 0./) | real |
| mid | Switch for mid-level convection<br>- 'NIL' : No mid-level convective scheme<br>- 'KF ' : Kain-Fritsch-based mid-level convective scheme | 'nil' | character(len=16) |
| mid_conserve | Conservation corrections for mid-level convective scheme<br>- 'NIL   ' : No conservation correction applied<br>- 'TEND  ' : Temperature and moisture tendencies corrected<br>- 'PRECIP' : Surface precipitation rate corrected | 'PRECIP' | character(len=16) |
| mid_depth | Minimum cloud depth for mid-level convection (m) | 2000. | real |
| mid_dpdd | Downdraft detrainment depth for mid-level convection (Pa) | 6000. | real |
| mid_emffrac | Fraction of environmental mass flux that enters updrafts | 'all' | character(len=16) |
| mid_emfmod | Modulation of the minimum environmental mass flux for mid-level convection | 'nil' | character(len=16) |
| mid_maxcape | Maximum deep CAPE (J/kg/m2) for mid-level convective triggering | -1 | real |
| mid_minbase | Minimum parcel departure level for mid-level convection (m) | 500. | real |
| mid_minemf | Minimum environmental mass flux for mid-level convection (kg/s) | 1e7 | real |
| mid_peff |  | -1. | real |
| shal | Switch for shallow convection<br>- 'NIL'<br>- 'KTRSNT'<br>- 'KTRSNT_MG'<br>- 'BECHTOLD' | 'nil' | character(len=16) |
| shal_conserve | Conservation corrections for shallow convective scheme<br>- 'NIL ' : No conservation correction applied<br>- 'TEND' : Temperature and moisture tendencies corrected | 'NIL' | character(len=16) |
| shal_timeconv |  | -1. | real |
| triglat | Over land and lakes we keep the value set by the "ramp" above over sea water:<br>- for :lat: >= TRIGLAT(2) we keep value set by the "ramp" KFCTRIG4<br>- for :lat: <= TRIGLAT(1) we use the new value KFCTRIGL<br>- and linear interpolation in between TRIGLAT(1) and TRIGLAT(2) | 0.0 | real |


### physics_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| acchr | Time length (hours) for special time accumulated physics variables | 0 | integer |
| advectke | Turbulent kinetic energy advect. is active if .true. | .false. | logical |
| chauf | Surface heat flux from oceans is active if .true. | .true. | logical |
| clip_tr_l | Clip tracers negative values | .true. | logical |
| cond_conserve | Conservation corrections for gridscale condensation<br>- 'NIL ' : No conservation correction applied<br>- 'TEND' : Temperature and moisture tendencies corrected | 'NIL' | character(len=16) |
| cond_evap | Evaporation parameter for Sunqvist gridscale condensation | 2.e-4 | real |
| cond_hmrst | Minimum cloud mixing ratio (kg/kg) for autoconversion in<br>Sunqvist gridscale condensation | 3.e-4 | real |
| cond_hu0max | Max allowed values of modified hu00 (threshold relative humidity<br>for stratiform condensation, Sunqvist gridscale condensation) | 0.975 | real |
| cond_hu0min | Min allowed values of modified hu00 (threshold relative humidity<br>for stratiform condensation, Sunqvist gridscale condensation) | 0.85 | real |
| cond_iceacc | Dry accretion of ice crystals by snow in Sundqvist (factor) | 5. | real |
| cond_infilter | Standard deviation length scale (gridpoints) of Gaussian smoother<br>applied to temperature and humidity inputs for Sunqvist gridscale<br>condensation) | -1. | real |
| debug_mem_l | Activate Debug memory mode | .false. | logical |
| debug_trace_l | Print a trace of the phy functions (MSG verbosity = debug) | .false. | logical |
| diag_twind | Obtain estimate of surface wind gusts if .true. | .false. | logical |
| diffuw | Diffuse vertical motion if .true. | .false. | logical |
| drag | Surface friction is active if .true.<br>Uses Schuman-Newell lapse rate if .false. | .true. | logical |
| etrmin2 | Minimal value for TKE in stable case (for 'CLEF') | 1.E-4 | real |
| evap | Surface evaporation is active if .true. | .true. | logical |
| fluvert | Boundary layer processes<br>- 'NIL    ': no vertical diffusion<br>- 'CLEF   ': non-cloudy boundary layer formulation<br>- 'MOISTKE': cloudy boundary layer formulation<br>- 'SURFACE': TODO<br>- 'SIMPLE ': a very simple mixing scheme for neutral PBLs | 'NIL' | character(len=16) |
| fnn_mask | (MOISTKE only) Apply factor fnn_reduc<br>- .false.: everywhere<br>- .true.: over water only | .false. | logical |
| fnn_reduc | (MOISTKE only) Reduction factor (between 0. and 1.) to be applied to the<br>parameter FNN (turbulent flux enhancement due to boundary layer clouds) | 1. | real |
| fnnmod | (CLEF+CONRES only) Non-dimensional parameter (must be >= 1.) that controls<br>the value of the flux enhancement factor in CONRES | 2. | real |
| fomic | Use Fomichev radiation code if .true. | .false. | logical |
| gwdrag | Gravity wave drag formulation<br>- 'NIL  ': no Gravity wave drag<br>- 'GWD86': gravity wave drag + low-level blocking<br>- 'SGO16': new formulation (2016) of GWD86 | 'NIL' | character(len=16) |
| hines_flux_filter | Number of times the 3-point filter will be applied to smooth the GW flux profiles | 0 | integer |
| iheatcal | Consider heating from non-orog. drag if = 1 | 0 | integer |
| indiag_list_s | Comma-separated list of diagnostic level inputs to read.<br>Default: indiag_list_s(1) = 'DEFAULT LIST',<br>expanded to: UU, VV, TT, HU + all dynamic Tracers | ' ' | character(len=32) |
| inilwc | Initialize water content and cloud fraction seen by radiation for time 0 if .true. | .false. | logical |
| input_type | Type of input system used<br>- 'DIST   ' : GEM 5.0 input system, RPN_COMM_IO/RPN_COMM_ezshuf_dist based<br>- 'BLOC   ' : GEM 5.0 input system, RPN_COMM_bloc based<br>- 'GEM_4.8' : GEM 4.8 input system, RPN_COMM_bloc based | 'DIST' | character(len=16) |
| intozot | Update ozone climatology during the run | .false. | logical |
| iuv_method | UV Method Optical properties of liquid cloud from condensation scheme for radiation<br>- 'NIL'        : No UV calculation<br>- 'IntegFit'   : Integration over the four UV bands with integrand weighting<br>- 'LinearFit'  : Fitted sum of the four UV broadband irradiances<br>- 'BandRatio'  : Scaling of input clear-sky UV index | 'IntegFit' | character(len=10) |
| kntrad_s | Time between full radiation calculation (units D,H,M,S,P) | '' | character(len=16) |
| kntraduv_s | Time between UV radiation calculation (units D,H,M,S,P) | '' | character(len=16) |
| kticefrac | Compute ice fraction in KTRSNT_MG if .true. | .true. | logical |
| lhn | use Latent Heat Nudging for the assmilation of radar-inferred precipitation rates<br>- 'NIL'   : No Latent Heat Nudging<br>- 'IRPCP' : Latent heat nudging from radar-inferred precipitation rates | 'NIL' | character(len=16) |
| lhn_filter | Standard deviation length scale (gridpoints) of Gaussian smoother<br>applied to RDPR, PR and TA before the application of Latent Heat Nudging<br>No smoothing applied when -ve | -1. | real |
| lhn_rainfix_l | Key to activate bugfix preventing inadvertent division of LHN tendencies by dt for precip-precip cases<br>bugfix will be active with:   lhn_rainfix_L = .true. | .false. | logical |
| lhn_weight | Modulation factor for the magnitude of Latent Heat Nudging being applied<br>modulated_tendencies = lhn_weight*(LHN tendencies) | 0. | real |
| linoz_chm | LINOZ prognostic stratospheric ozone<br>- 'NIL     ' :<br>- 'OZONE   ' :<br>- 'GHG     ' :<br>- 'OZONEGHG' : | 'NIL' | character(len=10) |
| lmetox | Add methane oxydation as source of humidity in the stratosphere if .true. | .false. | logical |
| longmel | Mixing length calc. scheme<br>- 'BLAC62  ': mixing length calc. using Blackadar<br>- 'BOUJO   ': mixing length calc. using Bougeault<br>- 'TURBOUJO': mixing length calc. using Bougeault in turbulent regimes (otherwise Blackadar)<br>- 'LH      ': mixing length calc. using Lenderink and Holtslag | 'BLAC62' | character(len=16) |
| moyhr | Time length (hours) for special time averaged physics variables | 0 | integer |
| mp_aeroact | Switch for aerosol activation scheme (1 = default, 2 = ARG + Aerosol climatology) | 1 | integer |
| my_ccntype | Switch for airmass type (1 = maritime, 2 = continental) | 1 | integer |
| my_dblmom_c | Double-moment for cloud (for 'mp_my' only) | .true. | logical |
| my_dblmom_g | Double-moment for graupel (for 'mp_my' only) | .true. | logical |
| my_dblmom_h | Double-moment for hail (for 'mp_my' only) | .true. | logical |
| my_dblmom_i | Double-moment for ice (for 'mp_my' only) | .true. | logical |
| my_dblmom_r | Double-moment for rain (for 'mp_my' only) | .true. | logical |
| my_dblmom_s | Double-moment for snow (for 'mp_my' only) | .true. | logical |
| my_diagon | Compute MY Diagnostic fields if .true. | .true. | logical |
| my_iceon | Ice-phase switched on if .true. | .true. | logical |
| my_initn | Initialize the number concentration for each category | .true. | logical |
| my_rainon | Autoconversion (cloud to rain) switched on | .true. | logical |
| my_sedion | Sedimentation switched on | .true. | logical |
| my_snowon | Snow initiation switched on | .true. | logical |
| my_tc3comp | Parameters for three-component freezing term | (/-5.,-5.,-5./) | real |
| my_warmon | Warm-phase switched on | .true. | logical |
| ninblocx | Physic input blocking along X | 1 | integer |
| ninblocy | Physic input blocking along Y | 1 | integer |
| non_oro | Hines non-orographic GWD scheme is active if .true. | .false. | logical |
| non_oro_pbot | Pressure (in Pa) that defines the bottom emission level for gravity waves | 61000.0 | real |
| nsloflux | Number of timesteps for which surface fluxes "FC" and "FV" are<br>gradually set from 0 to their full value in a "slow start fashion"<br>at the beginning of a time integration | 0 | integer |
| p3_clip_qv | switch for clipping of negative water vapor mixing ratios in microphysics (P3) | .true. | logical |
| p3_debug | switch for real-time debugging in microphysics (P3) | .false. | logical |
| p3_depfact | calibration factor for ice deposition in microphysics (P3) | 1.0 | real |
| p3_dtmax | Maximum time step (s) to be taken by the microphysics (P3) scheme, with time-splitting<br>used to reduce step to below this value if necessary | 60. | real |
| p3_ncat | Number of ice-phase hydrometeor categories to use in the P3 microphysics<br>scheme (currently limited to <5) | 1 | integer |
| p3_pfrac | precipitation fraction factor used by SCPF in microphysics (P3) | 1.0 | real |
| p3_resfact | model resolution factor used by SCPF in microphysics (P3) | 1.0 | real |
| p3_scpf_on | switch for subgrid cloud/precipitation fraction scheme (SCPF) in microphysics (P3) | .false. | logical |
| p3_subfact | calibration factor for ice sublimation in microphysics (P3) | 1.0 | real |
| p_runlgt | Vectoc lenght physics memory space folding for openMP | -1 | integer |
| pbl_cmu_timeavg | Time-averaging of transfer coefficient for momentum to reduce 2-dt<br>oscillations in fluxes | .false. | logical |
| pbl_conserve | Conservation corrections for PBL scheme<br>- 'NIL ' : No conservation correction applied<br>- 'TEND' : Temperature and moisture tendencies corrected | 'NIL' | character(len=16) |
| pbl_cucloud | Include the turbulent effects of trade wind cumulus clouds | .true. | logical |
| pbl_diff_condens | Diffuse condensate fields | .false. | logical |
| pbl_diss | Run with a modified closure for the dissipation length scale<br>- 'NIL  ' : No modified closure for the dissipation length scale<br>- 'LIM50' : A maximum value of 50m is imposed on dissipation length | 'NIL' | character(len=16) |
| pbl_dissheat | Dissipative heating tendencies are computed for the PBL scheme such<br>that total energy (kinetic + internal) is conserved<br>- 'NIL       ' : No dissipative heating is computed<br>- 'LOCAL_K   ' : Local total energy conservation based on diffusion coefficients<br>- 'LOCAL_TEND' : Local total energy conservation based on wind tendencies | 'NIL' | character(len=16) |
| pbl_flux_consistency | Call surface scheme immediately before the PBL | .false. | logical |
| pbl_func_stab | Class of stability functions (stable case) to use in the PBL<br>- 'DELAGE97  ' : Use functions described by Delage (1997; BLM)<br>- 'BELJAARS91' : Use functions described by Beljaars and Holtslag (1991; JAM)<br>- 'LOCK07    ' : Use functions described by Lock (2007; Tech Report) employed at UKMO | 'DELAGE97' | character(len=16) |
| pbl_func_unstab | Class of stability functions (unstable case) to use in the PBL<br>- 'DELAGE92' : Use functions described by Delage and Girard (1992; BLM)<br>- 'DYER74  ' : Use functions described by Dyer (1974; BLM) | 'DELAGE92' | character(len=16) |
| pbl_mlblac_max | Choose form of asymptotic mixing length for Blacadar-type estimates<br>- 'BLAC62' : Asymptotic 200 m proposed by Blackadar (1962; JGR) with clipping<br>- 'LOCK07' : Diagnosed asymptotic scale of Lock (2007; Tech Report) used at UKMO | 'BLAC62' | character(len=16) |
| pbl_mlturb_diss | Apply "turboujo" turbulence conditions to dissipation length scale | .false. | logical |
| pbl_moistke_legacy_cloud | Run with legacy moistke clouds (no limits on cloud effects) | .false. | logical |
| pbl_nonloc | Use the non-local PBL cloud formulation<br>- 'NIL   ' : no non-local PBL cloud formulation<br>- 'LOCK06' : Non-local cloud scheme of Lock and Mailhot (2006) | 'NIL' | character(len=16) |
| pbl_ribkg | Use the mixing length to average the Richardson number profile of (potentially)<br>many layers to derive a "background" Ri estimate | .false. | logical |
| pbl_ricrit | Richardson num. critical values for hysteresis | 1. | real |
| pbl_shal | PBL representation of boundary layer clouds<br>- 'NIL     ': No Shallow convection<br>- 'CONRES  ': Bulk Richardson number-based turbulent enhancement<br>- 'SHALOW  ': Deprecated (see 1998 RPN physics doc)<br>- 'SHALODQC': Deprecated (see 1998 RPN physics doc)<br>- 'GELEYN  ': Deprecated (see 1998 RPN physics doc) | 'NIL' | character(len=16) |
| pbl_slblend_layer | Layer over which to adjust from SL to PBL stability functions [(bot,top) in m] |  |  |
| pbl_tkediff | Adjustment to coefficient for TKE diffusion | 1. | real |
| pbl_tkediff2dt | Control of time scale for TKE diffusion | .false. | logical |
| pbl_turbsl_depth | Depth (Pa) of the always-turbulent near-surface layer in the PBL | 3000. | real |
| pbl_zerobc | Use true (motionless) surface boundary conditions for TKE diffusion | .false. | logical |
| pbl_zntau | Relaxation timescale (s) for mixing length smoothing | 7200. | real |
| pcptype | Scheme to determine precipitation type<br>- 'NIL     ': no call to bourge<br>- 'BOURGE  ': use Bourgouin algorithm (bourge1) to determine precip. types.<br>- 'BOURGE3D': | 'NIL' | character(len=16) |
| phystat_2d_l | Physic statistics output for 3d varables:<br>- .false. : mean, var, min and max for the whole 3d fiels<br>- .true.  : mean, var, min and max are done for each levels independently | .false. | logical |
| phystat_dble_l | Use double presision for physic statistics output | .false. | logical |
| phystat_freq_s | Physic statistics output Frequency | '0h' | character(len=16) |
| phystat_input_l | Print stats for phy_input read var | .false. | logical |
| phystat_list_s | Physic statistics output: bus variable list that should be included in physics<br>"block" stats. Possible values:<br>- Long varnames<br>- Short varnames<br>- 'ALLVARS=EDPV': all variables from E, D, P, V buses (any combination of the 4 letters); | ' ' | character(len=32) |
| qcfc11 | CFC11 bckgrnd atmospheric concentration (PPMV) | -1. | real |
| qcfc12 | CFC12 bckgrnd atmospheric concentration (PPMV) | -1 | real |
| qch4 | CH4 bckgrnd atmospheric concentration (PPMV) | -1. | real |
| qco2 | CO2 bckgrnd atmospheric concentration (PPMV) | -1. | real |
| qn2o | N2O bckgrnd atmospheric concentration (PPMV) | -1. | real |
| rad_atmpath | Atmospheric path length for solar radiation<br>- 'RODGERS67' : Formulation used by Li and Barker (2005)<br>- 'LI06' : Estimate of  Li and Shibata (2006) | 'RODGERS67' | character(len=16) |
| rad_cond_rei |  | -1. | real |
| rad_cond_rew |  | -1. | real |
| rad_conserve | Conservation corrections for radiation scheme<br>- 'NIL ' : No conservation correction applied<br>- 'TEND' : Temperature and moisture tendencies corrected | 'NIL' | character(len=16) |
| rad_esfc | Use emissivity computed by the surface schemes | .false. | logical |
| rad_linoz_l | Use LINOZ prognostic Ozone in radiation (CCCMARAD2 .and. LINOZ only) | .false. | logical |
| radfiles | format of radiation files to be read<br>- 'STD': RPN standard file<br>- 'UNF': unformatted | 'STD' | character(len=16) |
| radfix | Radiation fixes near the model top(for newrad only) if .true. | .true. | logical |
| radfltr | Vertical smoothing on radiative fluxes(for newrad only) if .true. | .true. | logical |
| radghg_l | Use climatological values of GHG in radiation (CCCMARAD2 only) | .false. | logical |
| radia | Radiation scheme<br>- 'NIL      ': no radiation scheme<br>- 'NEWRAD   ': complete radiation scheme<br>- 'CCCMARAD ': most advanced radiation scheme<br>- 'CCCMARAD2': most advanced radiation scheme v2 | 'NIL' | character(len=16) |
| radnivl | List of levels on which IR and VIS radiation calculations are<br>performed (to save on CPU time) (for newrad only) | 0 | integer |
| radslope | Key for activation of the radiation along slopes | .false. | logical |
| rmscon | Launching level value of GW RMS wind (m/s) from non-orographic origin | 1.0 | real |
| satuco | water/ice phase for saturation calc. if .true.;<br>water phase only for saturation calc. if .false. | .true. | logical |
| sfcflx_filter_iter | Surface fluxes figital filter, number of iterations | 1 | integer |
| sfcflx_filter_order | Surface fluxes digital filter order<br>- -1 : No filter<br>- 2 or 4: Apply filter 2nd or 4th order respectively | -1 | integer |
| sgo_bhfac | Tuning factor for blocking height | 1.5 | real |
| sgo_cdmin | Sets the minimum value of the drag coefficient in the orographic<br>blocking scheme. | 1.0 | real |
| sgo_nldirfac | Turns on/off the non-linear amplification factor (depending on wind<br>direction) of the drag coefficient in the orographic blocking scheme | .true. | logical |
| sgo_phic | Critical phase for blocking height | 0.2 | real |
| sgo_stabfac | Turns on/off the amplification factor (due to stability) of the drag<br>coefficient in the orographic blocking scheme | .true. | logical |
| sgo_tdfilter | Standard deviation length scale (gridpoints) of Gaussian smoother<br>applied to wind GWD tendencies | 1. | real |
| sgo_windfac | Description of threshold for mean wind speed for blocking |  |  |
| simisccp | (DEPRECATED) Run ISCCP cloud simulator (cccmarad only) if .true.<br>WARNING: This option is no longuer suppored, will be removed | .false. | logical |
| stcond | Condensation scheme name<br>- 'NIL       ' : No explicit condensation scheme used<br>- 'CONSUN    ' : Sunqvist type condensation scheme<br>- 'NEWSUND   ' : Sunqvist type condensation scheme<br>- 'MP_MY2_OLD' : Milbrandtl and Yau microphysics scheme (old formulation)<br>- 'MP_MY2    ' : Milbrandtl and Yau microphysics scheme<br>- 'MP_P3     ' : P3 microphysics scheme | 'NIL' | character(len=16) |
| stratos | Special treatment of stratosphere;<br>if .true. ignore convection/condensation tendencies where pressure is lower<br>than topc or specific humidity is lower than minq as specified in nocld.cdk | .false. | logical |
| taufac | Factor used in the gwd formulation = 1/(LENGTH SCALE) | 8.E-6 | real |
| test_phy | Run the physics in test harness mode | .false. | logical |
| timings_l | Print runtime timings | .false. | logical |
| tofd | Select a turbulent orographic form drag scheme<br>- 'NIL'        : No turbulent orographic form drag scheme<br>- 'BELJAARS04' : Form drag scheme described by Beljaars et al. (2006; QJRMS)<br>WARNING: This option is broken thus disabled- will be fixed in dev branch | 'NIL' | character(len=16) |
| tofd_alpha | BELJAARS04 turbulent orographic form drag scheme alpha parameter | 12. | real |
| ts_flxir | (newrad only) Use TT(12000) instead of skin temp in downward IR<br>flux calculation if .true. | .false. | logical |


### series Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| p_serg_serstp | Times series package stops at this timestep | huge(1) | integer |
| p_serg_srprf_s | List of time series for profile variables | ' ' | character(len=SER_STRLEN_VAR) |
| p_serg_srsrf_s | List of time series for surface variables | ' ' | character(len=SER_STRLEN_VAR) |
| p_serg_srwri | Number of timesteps between time-series writeout | 1 | integer |
| xst_stn_latlon | Stations chosen in lat,lon for time-series<br>Format: "STN1_NAME",lat1,lon1, "STN2_NAME",lat2,lon2, ... |  |  |


### surface_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| adj_i0_snow | Adjust surface temperature over snow after reading (coherency check) | .true. | logical |
| beta | Prandtl number for neutral stability (initialized by SL module) | 0. | real |
| diusst | Diurnal SST scheme<br>- 'NIL    ' : No Diurnal SST scheme<br>- 'FAIRALL' : #TODO: define | 'NIL' | character(len=16) |
| diusst_coolskin | Diurnal SST scheme active coolskin if .true. | .true. | logical |
| diusst_coolskin_lakes | Diurnal SST scheme active coolskin over freshwater lakes if .true. | .true. | logical |
| diusst_warmlayer | Diurnal SST scheme active warmlayer if .true. | .true. | logical |
| diusst_warmlayer_lakes | Diurnal SST scheme active warmlayer over freshwater lakes if .true. | .true. | logical |
| dp_svs | Depth of soil layers in [METERS] in SVS land surface scheme (schmsol=SVS) | -1.0 | real |
| ice_emiss |  | -1. | real |
| icelac | Set water temperature of ice-covered lakes to 0C for points north of<br>ice line if .true.<br>needs an initialization file otherwise the model stops | .false. | logical |
| icemelt | Sea ice melting | .false. | logical |
| impflx | Implicit surface fluxes if .true.; explicit fluxes if .false. | .false. | logical |
| isba_melting_fix | If .true. apply temporary fix to ISBA<br>- timestep dependent KCOEF<br>- No PSN factor for meting and freezing | .false. | logical |
| isba_snow_z0veg | Use the vegetation-only roughness length to compute vegetation snow fraction | .false. | logical |
| isba_soil_emiss |  | -1. | real |
| isba_zr_freeze | If .true., freeze precipitation reaching the ground in sub-zero conditions | .false. | logical |
| kdp | OBSOLETE, REPLACED by KHYD !!! WILL BE EVENTUALLY REMOVED<br>Deepest active (permeable) soil layer in SVS land surface scheme (schmsol=SVS) | -1 | integer |
| khyd | Last/Deepest soil layer considered during the accumulation of<br>lateral flow and drainage. Drainage is taken as the vertical flux<br>leaving layer KHYD, and lateral flow as the sum of lateral flows from<br>layers 1 to KHYD | -1 | integer |
| kntveg_s | Vegetation field update frequency (units D,H,M,S,P) | '' | character(len=16) |
| lake_leadfrac | Lead fraction for ice-covered lakes | 0. | real |
| leadfrac | Minimum fraction of leads in sea ice.&nbsp; Multiply ice fraction by (1.-leadfrac) | 0.03 | real |
| limsnodp | Limit snow depth to 10 cm for calculation of heat conductivity of snow<br>over sea-ice and glacier if .true. | .false. | logical |
| owflux | (coupling) fluxes over ocean are taken from ocean model if .true. | .false. | logical |
| read_emis | read-in land surface emissivity if .true. | .false. | logical |
| read_z0vh | read-in high vegetation roughness for SVS if .true. | .false. | logical |
| salty_qsat | Takes into account effect of ocean salinity on saturation specific<br>humidity at ocean surface (boundary condition for LH flux calculation) | .false. | logical |
| schmsol | Land surface processes<br>- 'NIL ' : No Land surface processes<br>- 'ISBA' : Interaction Soil Biosphere Atmosphere (ISBA) land sfc scheme<br>- 'SVS ' : Soil, Vegetation, and Snow (SVS) (Multibudget) land sfc scheme | 'ISBA' | character(len=16) |
| schmurb | Urban surface processes<br>- 'NIL' : No Urban surface processes<br>- 'TEB' : Town Energy Balance (TEB) urban scheme | 'NIL' | character(len=16) |
| sl_func_stab | Class of stability functions (stable case) to use in the surface layer<br>- 'DELAGE97  ' : Use functions described by Delage (1997; BLM)<br>- 'BELJAARS91' : Use functions described by Beljaars and Holtslag (1991; JAM)<br>- 'LOCK07    ' : Use functions described by Lock (2007; Tech Report) employed at UKMO | 'DELAGE97' | character(len=16) |
| sl_func_unstab | Class of stability functions (unstable case) to use in the surface layer<br>- 'DELAGE92' : Use functions described by Delage and Girard (1992; BLM)<br>- 'DYER74  ' : Use functions described by Dyer (1974; BLM) | 'DELAGE92' | character(len=16) |
| sl_lmin_glacier | Minimum Obukhov length (L) for glaciers | -1. | real |
| sl_lmin_seaice | Minimum Obukhov length (L) for sea ice | -1. | real |
| sl_lmin_soil | Mimimum Obukhov length (L) for soil surfaces | -1. | real |
| sl_lmin_water | Minimum Obukhov length (L) for water | -1. | real |
| sl_rineutral | Define bulk Ri values for near-neutral regime in the surface layer | 0. | real |
| sl_z0ref | Use a reference roughness for surface layer calculations | .false. | logical |
| snoalb_anl | Use snow albedo "I6" directly if .true.;<br>Use snow age "XA" to calculate snow albedo if .false. | .true. | logical |
| snow_emiss |  | -1. | real |
| soiltext | Soil texture database/calculations for SVS land surface scheme<br>- 'GSDE   '   : 8 layers of sand & clay info from Global Soil Dataset for ESMs (GSDE)<br>- 'SLC    '   : 5 layers of sand & clay info from Soil Landscape of Canada (SLC)<br>- 'SOILGRIDS' : 7 layers of sand & clay info from ISRIC ? World Soil Information | 'GSDE' | character(len=16) |
| svs_dynamic_z0h | use dynamic calculation of z0h for bare ground + vegetation  for SVS if .true. | .false. | logical |
| svs_hrsurf_sltext | use hrsurf based on soil texture for SVS if .true. | .false. | logical |
| svs_local_z0m | use local momentum (no snow) roughness for SVS if .true. | .false. | logical |
| tdiaglim | Limit temperature inversions to 8K/40m in surface layer if .true. | .false. | logical |
| urban_params_new | New urban surface parameters | .false. | logical |
| use_eff_surf_tq |  | .false. | logical |
| use_photo | OPTION TO USE PHOTOSYNTHESIS CODE FOR STOMATAL RESISTANCE in SVS | .true. | logical |
| veg_rs_mult | Factor multiplying stomatal resistance in ISBA | 1. | real |
| vf_type | VF definitions and mapping in SVS<br>- 'CLASSIC' : Same VF definitions as ISBA<br>- 'CCILCECO' : New VF definitions used in SVS only<br>with geo. fields generated using CCILC 2015 + Ecobiomes | 'CLASSIC' | character(len=16) |
| water_emiss |  | -1. | real |
| z0dir | Use directional roughness length if .true. | .false. | logical |
| z0hcon | Constant value of thermal roughness length (m) applied over water within<br>latitudinal band defined by z0tlat | 4.0e-5 | real |
| z0min | Minimum value of momentum roughness length (m) | 1.5e-5 | real |
| z0mtype | Momentum roughness length formulation over water<br>- 'CHARNOCK' : #TODO: define<br>- 'BELJAARS' : #TODO: define | 'CHARNOCK' | character(len=16) |
| z0seaice | Roughness length for sea ice | 1.6e-4 | real |
| z0tevol | Thermal roughness length formulation over vegetation<br>- 'FIXED' : Uses z0h = z0m<br>- 'ZILI95': evolves with u* | 'FIXED' | character(len=16) |
| z0tlat | Latitude (2 elements, in degrees) used to specify Z0T over water<br>- If :lat: <= Z0TLAT(1) constant Z0T.<br>- If :lat: >= Z0TLAT(2) Charnock's relation.<br>- In between, linear interpolation is used. | 0. | real |
| z0ttype | Thermal roughness length formulation over water<br>- 'MOMENTUM' : Uses z0h = z0m (replaces key z0trdps300=.false.)<br>- 'DEACU12'  : #TODO: define  (replaces key z0trdps300=.true.)<br>- 'ECMWF'    : #TODO: define  (New formulation used by ECMWF) | 'MOMENTUM' | character(len=16) |
| zt | Height at which to compute screen-level temperature (m) | 1.5 | real |
| zu | Height at which to compute anemomenter-level winds (m) | 10. | real |


