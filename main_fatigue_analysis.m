%% MAIN_FATIGUE_ANALYSIS  Time-domain random-wave fatigue of monopile welds.
% =========================================================================
%  ACCRAR offshore wind farm - life-extension study (Imperial College GDP /
%  Rambull). Quasi-static, WAVE-ONLY cumulative fatigue damage at Nw = 37 welds,
%  following the workflow: JONSWAP -> random-phase eta(t) -> Airy kinematics
%  (vertical extrapolation) -> Morison force -> per-weld bending moment ->
%  hot-spot stress -> built-in rainflow -> Miner with a bilinear S-N curve.
%
%  Files in this package:
%    main_fatigue_analysis.m        (this driver)
%    generate_jonswap_timeseries.m
%    calculate_morison_moment.m     (now returns M for an ARRAY of welds)
%    calculate_damage_rainflow.m
%
%  Requires: MATLAB R2016b+ and RAINFLOW (Signal Processing Toolbox, R2017b+).
%
%  EFFICIENCY (Requirement 3): the wave kinematics over the water column are
%  IDENTICAL for every weld in a given sea state, so they are generated ONCE per
%  sea state. calculate_morison_moment then returns M as an (Nw x Nt) matrix --
%  one moment history per weld -- by integrating the shared force field from each
%  weld's elevation up to eta. Only rainflow/Miner (which cannot be vectorised
%  across signals) is looped per weld.
%
%  READ BEFORE TRUSTING ANY NUMBER (see notes at the bottom of this file):
%    1. The 37 WELD properties (z_weld, OD_steel, t_wall, SCF) are PLACEHOLDERS.
%       Overwrite the four arrays below with the real per-weld values (drawing
%       can table + confirmed SCFs) before reporting.
%    2. The WAVE SCATTER is CONSTRUCTED, not measured: the Design Basis contains
%       no long-term/operational scatter, only extreme return-period criteria.
%       See the heavily-commented scatter_data block for exactly what was used.
%    3. This is QUASI-STATIC and WAVE-ONLY. Real monopile fatigue is driven by
%       structural DYNAMICS (1st mode ~0.2 Hz, between 1P=0.129 and 3P=0.387 Hz)
%       and by WIND (Table 9-3 FDEL). SDOF amplification is provided but OFF.
% =========================================================================

clear; clc;

%% ----------------------------- SITE / PHYSICS ---------------------------
d        = 63.2;          % mean water depth, MSL -> mudline      [m]
g        = 9.81;          % gravity                                [m/s^2]
rho      = 1025;          % seawater density                       [kg/m^3]
gamma_js = 3.3;           % JONSWAP peak-enhancement (North Sea)   [-]

%% --------------------------- SIMULATION CONTROL -------------------------
dur     = 3600;           % duration per sea-state realisation     [s] (1 hr)
dt      = 0.5;            % time step (Nyquist 1 Hz >> wave band)   [s]
Nfreq   = 2048;           % JONSWAP frequency components            [-]
Nz      = 80;            % vertical integration nodes (shared grid)[-]
Nseeds  = 1;             % random-phase seeds per sea state        [-] (raise for variance)

%% ----------------------------- GEOMETRY ---------------------------------
% Steel outer diameter vs elevation [mMSL]. Submerged shaft 11.0 m; tapers to
% 8.0 m near the pile top (~+5 mMSL). Clamped outside the defined range.
z_od   = [-63.2, -32.381, -7.5,2.4,  5.0];
OD_od  = [ 11.0, 11, 8, 8.0,  8.0];
geom.ODfun = @(z) interp1(z_od, OD_od, min(max(z,-91.6),5.0), 'linear');

% Marine-growth thickness vs elevation [mMSL] (UPDATED Table 12-3):
%   150 mm for +2 .. -40 mMSL, 75 mm below -40 mMSL, none above +2.
geom.mgfun = @(z) 0.150*(z>=-40 & z<=2) + 0.075*(z<-40);

%% ===================== WELDS: 37 LOCATIONS (PLACEHOLDERS) ===============
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Replace ALL FOUR arrays with the real per-weld values once available.
% Each is a column vector of length Nw = 37 (one entry per girth weld).
Nw       = 28;
z_weld   = [18.7
15.5
12
8.5
5.3
4.85
1
-3
4.76
2.5
-0.5
-4
-7.5
-11.06
-14.62
-18.18
-21.74
-25.3
-28.86
-32.381
-36.261
-40.141
-44.021
-47.901
-51.781
-55.661
-59.541
-63.2];     % weld elevations [mMSL], mudline -> top
OD_steel = geom.ODfun(z_weld);           % steel OD at each weld [m]      (placeholder)
t_wall   = [0.11
0.11
0.11
0.11
0.11
0.06
0.06
0.06
0.14
0.11
0.08
0.08
0.09
0.09
0.07
0.07
0.08
0.08
0.07
0.07
0.07
0.08
0.08
0.08
0.08
0.08
0.08
0.08];             % wall thickness at each weld [m] (placeholder)
SCF      = [1.00
1.07
1.13
1.06
1.00
1.00
1.00
1.00
1.00
1.14
1.17
1.02
1.30
1.00
1.18
1.10
1.00
1.04
1.02
1.43
1.02
1.04
1.06
1.04
1.02
1.17
1.23
1.02];             % total SCF at each weld [-]     (placeholder)
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% --------------------------- HYDRODYNAMICS ------------------------------
hydro.rho     = rho;
hydro.Cd      = 0.952;     % drag coefficient   [-]  (rough cylinder; tune per RP-C205)
hydro.Cm      = 2.00;     % inertia coefficient[-]  (rough cylinder; tune per RP-C205)
hydro.current = @(z) zeros(size(z));   % current OFF for the fatigue baseline.
% To include current, use the 1/7th power-law profile from Table 10-4, e.g.:
Us = 0.79;
hydro.current = @(z) Us*max((z+d)/d,0).^(1/7);

%% --------------------- SECTION (common) & S-N CURVE ---------------------
sec.t_ref = 0.025;        % reference thickness                  [m]
sec.k_thk = 0.05;         % thickness exponent                   [-]
% (per-weld OD_steel, t_wall, SCF are injected inside the loop)

% Bilinear DNV-RP-C203 F-curve, seawater with cathodic protection.
sn.m1    = 3.5;      sn.loga1 = 13.473;     % high-stress branch
sn.m2    = 5;      sn.loga2 = 16.377;     % low-stress branch
sn.Nknee = 5e6;                            % cycles at slope change (knee ~65.8 MPa)

DFF      = 3.0;           % design fatigue factor (Sec 6.4)        [-]
years_op = 25;            % operational design life (Sec 6 / Table 5-1) [yr]

%% ---- DYNAMIC AMPLIFICATION (optional, OFF by default) ------------------
opt.Nz       = Nz;
opt.kin      = 'extrapolated'; % splash-zone kinematics: DNV-RP-C205 [3.3.3.4] for linear seas
opt.dt       = dt;
opt.applyDAF = false;     % set true to switch on SDOF amplification
opt.fn       = 0.22;      % !!! PLACEHOLDER 1st natural frequency [Hz] - use PLAXIS NFA !!!
opt.zeta     = 0.04;      % !!! PLACEHOLDER structural damping ratio [-] !!!

%% ===================== WAVE FATIGUE SCATTER (CONSTRUCTED) ===============
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% The Design Basis (Design_Basis_GDP_2026_v1.0.pdf) does NOT contain a long-term
% / operational wave scatter diagram. Section 10.2 "Wave Conditions" (Table 10-1,
% p.13) and Section 12.1 (updated Table 12-2, p.16) give ONLY extreme omni-
% directional return-period criteria. The matrix below is therefore a CONSTRUCTED,
% representative fatigue scatter, built as follows:
%
%   * Hs-Tp relation: fitted to the DB extreme (Hs,Tp) pairs in Tables 10-1 and
%     12-2. Those pairs imply a near-constant peak steepness s_p = 2*pi*Hs/(g*Tp^2)
%     ~ 0.037-0.044, i.e.  Tp ~= 4.0*sqrt(Hs)  (s_p ~ 0.04). Tp values below are
%     taken from this relation. [DB Tables 10-1, 12-2]
%   * Upper Hs bound: anchored to the 1-yr return Hs = 6.3 m (updated Table 12-2,
%     Section 12.1). Operational fatigue is dominated by frequent low-moderate
%     seas, so bins run 0.5 -> 7.0 m (7.0 just above the 1-yr value, rare).
%   * The 50-yr / 100-yr seas (Hs = 11.0 / 13.5 m, Table 12-2) are ULS, NOT FLS,
%     and their true occurrence fraction is ~1e-5..1e-6; they are deliberately
%     EXCLUDED from this fatigue scatter.
%   * Probabilities (col 3): an ASSUMED monotonic-decaying marginal typical of an
%     energetic east-Scotland site, normalised to sum to 1.0. These occurrence
%     weights are NOT from the DB (it has no occurrence data) and are the single
%     biggest assumption here -- REPLACE with the site lumped scatter when available.
%
%   Columns:   Hs [m]   Tp [s] (=4.0*sqrt(Hs))   Probability [-]
%% ===================== WAVE FATIGUE SCATTER (JONSWAP METHOD) ============

db_extremes = [6.3, 10.4; 11.0, 12.9; 13.5, 14.0];   % DB Table 12-2

scatter_data = build_scatter_jonswap(g, gamma_js, db_extremes);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Hs_v   = scatter_data(:,1);
Tp_v   = scatter_data(:,2);
prob_v = scatter_data(:,3);
prob_v = prob_v / sum(prob_v);              % enforce normalisation
Nstate = size(scatter_data,1);

secPerYear = 365.25*24*3600;                % seconds per year [s]
reps_per_state = prob_v * secPerYear / dur; % how many 'dur'-records per year per state

%% ============================== MAIN LOOP ===============================
% Outer loop: sea states. Kinematics generated ONCE per (state, seed); the
% 37 weld moments come back together as M (Nw x Nt). Inner loop: per-weld
% rainflow/Miner damage. Damage accumulates into D_total (Nw x 1), in 1/yr.
D_total = zeros(Nw,1);

fprintf('\n  Running %d sea states x %d seed(s), %d welds, %.0f s @ dt=%.2g s ...\n', ...
        Nstate, Nseeds, Nw, dur, dt);
fprintf('  state |  Hs   Tp  | prob   |  reps/yr  | max sigma over welds [MPa]\n');
fprintf('  ------+-----------+--------+-----------+---------------------------\n');

for is = 1:Nstate
    Hs = Hs_v(is);  Tp = Tp_v(is);

    D_seed_w = zeros(Nw, Nseeds);           % per-weld damage per record, each seed
    smax_rec = 0;
    for js = 1:Nseeds
        seed = 1000*is + js;                % reproducible, distinct per (state,seed)

        % 1) random-phase JONSWAP surface elevation + components  (ONCE)
        [t, eta, comp] = generate_jonswap_timeseries(Hs, Tp, gamma_js, ...
                                                      dur, dt, Nfreq, d, seed);
        % 2) Airy kinematics (extrapolated) + Morison -> ALL 37 weld moments (ONCE)
        M = calculate_morison_moment(comp, t, eta, geom, hydro, z_weld, d, opt);  % Nw x Nt

        % 3) per-weld hot-spot stress + rainflow + Miner damage
        for i = 1:Nw
            sec.OD_steel = OD_steel(i);
            sec.t_wall   = t_wall(i);
            sec.SCF      = SCF(i);
            [Drec, info] = calculate_damage_rainflow(M(i,:), sec, sn, dt);
            D_seed_w(i,js) = Drec;
            smax_rec = max(smax_rec, info.sigma_max);
        end
    end

    D_state_w = mean(D_seed_w, 2);          % Nw x 1  mean damage per record over seeds
    D_total   = D_total + D_state_w * reps_per_state(is);   % accumulate annual damage

    fprintf('   %4d | %4.1f %5.1f | %5.3f  | %9.2e | %8.2f\n', ...
            is, Hs, Tp, prob_v(is), reps_per_state(is), smax_rec);
end

%% ============================== RESULTS =================================
% D_total is the annual damage per weld [1/yr].
life_yr  = 1 ./ D_total;                    % unfactored fatigue life per weld [yr]
util     = D_total * years_op * DFF;        % factored utilisation over operational life
pass_w   = util <= 1.0;
[util_max, igov] = max(util);               % governing (most-damaged) weld

fprintf('\n  =============================== PER-WELD FATIGUE SUMMARY ===============================\n');
fprintf('  weld |  z_weld | OD_st | t_wall |  SCF  |  D_annual   |  life [yr]  | util(=D*%2dy*%.1f) | check\n', years_op, DFF);
fprintf('  -----+---------+-------+--------+-------+-------------+-------------+------------------+------\n');
for i = 1:Nw
    if pass_w(i), chk = 'PASS'; else, chk = 'FAIL'; end
    fprintf('   %3d | %7.2f | %5.2f | %6.4f | %5.2f | %11.3e | %11.1f | %16.3e | %s\n', ...
            i, z_weld(i), OD_steel(i), t_wall(i), SCF(i), D_total(i), life_yr(i), util(i), chk);
end
fprintf('  =======================================================================================\n');

fprintf('\n  Governing weld: #%d at z = %.2f mMSL\n', igov, z_weld(igov));
fprintf('    annual damage ....... %12.4e  [1/yr]\n', D_total(igov));
fprintf('    unfactored life ..... %12.1f  [yr]\n', life_yr(igov));
fprintf('    utilisation (DFF=%.1f) %12.4e  ->  %s (<=1 required)\n', DFF, util_max, ternary(util_max<=1,'PASS','FAIL'));
fprintf('    welds failing ....... %d of %d\n', sum(~pass_w), Nw);

fprintf(['\n  REMINDERS: (1) weld arrays (z/OD/t/SCF) and the wave scatter are placeholders/constructed\n', ...
         '                 - replace with real per-weld data and the site scatter before reporting.\n', ...
         '             (2) wave-only & quasi-static - add structural dynamics (opt.applyDAF) and the\n', ...
         '                 wind FDEL (Table 9-3) via Miner combination for a defensible life number.\n\n']);

% =========================================================================
%  ENGINEER'S NOTES
% -------------------------------------------------------------------------
%  * Governing location: the largest wave bending moment is at the mudline, so the
%    lowest welds dominate D_total; welds high on the shaft see almost no wave
%    moment (their fatigue is wind/dynamics driven). This is the expected pattern.
%  * Efficiency: kinematics are O(Nz*Nfreq*Nt) and are computed ONCE per sea state
%    inside calculate_morison_moment; the 37 weld moments are 37 cheap masked
%    trapz integrations of the shared force field. Cost ~ independent of Nw.
%  * Dynamics (the big one): set opt.applyDAF = true and plug opt.fn from the
%    PLAXIS / NFA result and a measured/assumed opt.zeta. Quasi-static under-
%    predicts because the 1st mode (~0.2 Hz) sits inside the wave band.
%  * Large diameter: with D up to 11 m, short components approach the diffraction
%    regime (pi D/L -> 0.5+) where a constant Cm over-predicts inertia; the rigorous
%    fix is the MacCamy-Fuchs effective inertia coefficient for the HF tail.
%  * Convergence: re-run with larger Nfreq, smaller dt, more seeds and higher Nz
%    and confirm D_total (especially at the governing weld) is stable before quoting.
% =========================================================================

% --- tiny inline ternary helper -----------------------------------------
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
