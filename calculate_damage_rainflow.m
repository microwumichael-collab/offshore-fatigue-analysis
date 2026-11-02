function [D, info] = calculate_damage_rainflow(M, sec, sn, dt)
%CALCULATE_DAMAGE_RAINFLOW  Hot-spot stress, rainflow counting and Miner damage.
%
%   [D, info] = CALCULATE_DAMAGE_RAINFLOW(M, sec, sn, dt)
%
%   Converts a bending-moment time series at the weld into a hot-spot stress
%   series, rainflow-counts the stress ranges with MATLAB's built-in RAINFLOW,
%   and accumulates Palmgren-Miner damage against a bilinear S-N curve.
%
%   INPUTS
%     M     bending moment at the weld          [N*m]  (1 x Nt)
%     sec   struct describing the cross-section / hot spot:
%             .OD_steel  steel outer diameter at weld   [m]
%             .t_wall    (corroded) wall thickness        [m]
%             .SCF       total stress concentration factor [-]
%             .t_ref     reference thickness               [m]   (0.025)
%             .k_thk     thickness exponent                [-]   (0.20)
%     sn    struct describing the bilinear S-N curve (log10 form):
%             .m1     slope, high-stress branch         [-]   (3)
%             .loga1  intercept, high-stress branch     [-]   (11.455)
%             .m2     slope, low-stress branch          [-]   (5)
%             .loga2  intercept, low-stress branch      [-]   (15.091)
%             .Nknee  cycles at the slope change        [-]   (1e6)
%     dt    time step of M (for cycle-rate reporting)   [s]
%
%   OUTPUTS
%     D     cumulative Miner damage for this record     [-]
%     info  struct with diagnostics:
%             .sigma_MPa   hot-spot stress series        [MPa] (1 x Nt)
%             .ranges      rainflow stress ranges        [MPa]
%             .counts      rainflow cycle counts         [-]
%             .W           section modulus               [m^3]
%             .dsig_knee   S-N knee stress range          [MPa]
%             .n_cycles    total counted cycles           [-]
%             .sigma_max   peak |stress| in the record    [MPa]
%
%   Note on the thickness factor: DNV applies a penalty only for t > t_ref, so
%   the factor uses max(t_wall, t_ref) and is therefore >= 1.
% -------------------------------------------------------------------------

    if exist('rainflow','file') ~= 2
        error(['MATLAB''s built-in RAINFLOW was not found. It ships with the ', ...
               'Signal Processing Toolbox (R2017b+). Install/enable it, or ', ...
               'substitute an equivalent rainflow counter.']);
    end

    % --- section modulus of the (corroded) steel tube --------------------
    Do = sec.OD_steel;                          % outer steel diameter [m]
    Di = Do - 2*sec.t_wall;                      % inner diameter        [m]
    W  = pi*(Do^4 - Di^4) / (32*Do);             % elastic section modulus [m^3]

    % --- thickness correction factor (>= 1) ------------------------------
    tcf = (max(sec.t_wall, sec.t_ref)/sec.t_ref)^sec.k_thk;

    % --- hot-spot stress series  sigma = (M/W)*SCF*tcf -------------------
    sigma_Pa  = (M(:)/W) * sec.SCF * tcf;        % Pa  (column)
    sigma_MPa = sigma_Pa / 1e6;                  % MPa

    % --- rainflow counting (built-in) ------------------------------------
    % rainflow returns columns: [Count, Range, Mean, Start, End].
    rf     = rainflow(sigma_MPa);
    counts = rf(:,1);                            % cycle counts (may be 0.5 = half cycles)
    ranges = rf(:,2);                            % stress ranges [MPa]

    % --- bilinear S-N: cycles to failure for each range ------------------
    N = snLife(ranges, sn);

    % --- Palmgren-Miner accumulation -------------------------------------
    D = sum(counts ./ N);

    % --- diagnostics ------------------------------------------------------
    info.sigma_MPa = sigma_MPa.';
    info.ranges    = ranges;
    info.counts    = counts;
    info.W         = W;
    info.dsig_knee = 10^((sn.loga1 - sn.loga2)/(sn.m1 - sn.m2));
    info.n_cycles  = sum(counts);
    info.sigma_max = max(abs(sigma_MPa));
    if nargin >= 4 && ~isempty(dt)
        info.cycles_per_s = info.n_cycles / (numel(M)*dt);
    end
end

% =========================================================================
function N = snLife(dsig, sn)
%SNLIFE  Cycles to failure on a bilinear (two-slope) S-N curve, log10 form.
%   High-stress branch (slope m1) above the knee; low-stress branch (slope m2)
%   below it. The knee range is the intersection of the two branches, which is
%   continuous at N = Nknee by construction of the DNV constants.
    dsig_knee = 10^((sn.loga1 - sn.loga2)/(sn.m1 - sn.m2));   % MPa
    N = zeros(size(dsig));

    hi = dsig >  dsig_knee;                      % steep branch (m1)
    lo = ~hi;                                     % shallow branch (m2)

    N(hi) = 10.^(sn.loga1 - sn.m1*log10(dsig(hi)));
    N(lo) = 10.^(sn.loga2 - sn.m2*log10(dsig(lo)));

    % zero-range cycles (flat segments) contribute no damage -> infinite life
    N(dsig <= 0) = Inf;
end
