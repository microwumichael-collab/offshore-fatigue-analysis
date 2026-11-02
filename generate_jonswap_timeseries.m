function [t, eta, comp] = generate_jonswap_timeseries(Hs, Tp, gamma, dur, dt, Nfreq, d, seed)
%GENERATE_JONSWAP_TIMESERIES  Random-phase JONSWAP sea-surface elevation.
%
%   [t, eta, comp] = GENERATE_JONSWAP_TIMESERIES(Hs, Tp, gamma, dur, dt, Nfreq, d, seed)
%
%   Builds a 1- or 3-hour irregular sea-surface elevation eta(t) by linear
%   superposition of N harmonic components whose amplitudes are drawn from a
%   JONSWAP spectrum (DNV-RP-C205) and whose phases are random. Also returns
%   the per-component data (w, k, a, phi) so the kinematics can be built later.
%
%   INPUTS
%     Hs     significant wave height          [m]
%     Tp     peak (spectral) period           [s]
%     gamma  JONSWAP peak-enhancement factor   [-]   (3.3 = North Sea standard)
%     dur    simulation duration              [s]
%     dt     time step                        [s]
%     Nfreq  number of frequency components   [-]
%     d      mean water depth (MSL -> mudline)[m]
%     seed   RNG seed (reproducibility)       [-]
%
%   OUTPUTS
%     t      time vector                      [s]   (1 x Nt)
%     eta    surface elevation at the pile    [m]   (1 x Nt)
%     comp   struct of component data:
%              .w   angular frequencies        [rad/s]  (Nfreq x 1)
%              .k   wave numbers               [rad/m]  (Nfreq x 1)
%              .a   component amplitudes       [m]      (Nfreq x 1)
%              .phi random phases              [rad]    (Nfreq x 1)
%
%   Method: DNV-RP-C205 Sec.3.5 (JONSWAP) + linear (Airy) dispersion.
%   Requires MATLAB R2016b+ (implicit expansion).
% -------------------------------------------------------------------------

    g = 9.81;
    wp = 2*pi/Tp;                         % peak angular frequency [rad/s]

    % --- FIXED absolute frequency band (independent of Tp) -----------------
    % Chosen wide enough to capture all JONSWAP energy for Tp in [5,14] s,
    % and fine enough that the record repeat period exceeds 'dur'.
    wmin = 0.10;                          % [rad/s]  ~ T = 63 s
    wmax = 2.50;                          % [rad/s]  ~ T = 2.5 s
    dw   = (wmax - wmin)/Nfreq;
    w    = (wmin + ((1:Nfreq)' - 0.5)*dw);% mid-bin angular frequencies (col)

    % --- repeat-period guard ----------------------------------------------
    Trepeat = 2*pi/dw;                    % record repeats every Trepeat seconds
    if Trepeat < dur
        warning('jonswap:repeat', ...
          ['Record repeat period (%.0f s) < duration (%.0f s). ', ...
           'Increase Nfreq to avoid spurious periodicity.'], Trepeat, dur);
    end

    % --- JONSWAP spectrum (DNV-RP-C205) -----------------------------------
    sigma          = 0.07*ones(size(w));  % spectral width parameter
    sigma(w > wp)  = 0.09;
    Ag  = 1 - 0.287*log(gamma);           % normalising factor A_gamma
    Spm = (5/16) * Hs^2 * wp^4 .* w.^(-5) .* exp(-1.25*(w/wp).^(-4));   % P-M part
    S   = Ag * Spm .* gamma.^( exp( -0.5*((w - wp)./(sigma.*wp)).^2 ) );% [m^2 s/rad]

    % --- deterministic-amplitude / random-phase components ----------------
    a = sqrt(2*S*dw);                     % component amplitudes [m]
    rng(seed);
    phi = 2*pi*rand(Nfreq,1);             % random phases [rad]

    % --- wave numbers from linear dispersion (vectorised) -----------------
    k = linearDispersion(w, d, g);

    % --- surface elevation time series ------------------------------------
    t     = (0:dt:dur);                   % 1 x Nt
    Theta = w*t + phi;                    % Nfreq x Nt
    eta   = (a.') * cos(Theta);           % 1 x Nt   ( sum_i a_i cos(theta_i) )

    % --- sanity: recovered Hs from spectral variance ----------------------
    m0     = sum(S)*dw;                   % zeroth spectral moment
    Hs_chk = 4*sqrt(m0);
    if abs(Hs_chk - Hs)/Hs > 0.02
        warning('jonswap:variance', ...
          'Recovered Hs = %.2f m vs target %.2f m (widen band / raise Nfreq).', Hs_chk, Hs);
    end

    comp = struct('w', w, 'k', k, 'a', a, 'phi', phi);
end

% =========================================================================
function k = linearDispersion(w, d, g)
%LINEARDISPERSION  Solve w^2 = g k tanh(k d) for k (vectorised, Newton).
    p = w.^2 * d / g;                     % Fenton & McKee initial guess
    q = tanh(p.^0.75).^(-2/3);
    k = w.^2 .* q / g;
    for it = 1:8
        f  = g*k.*tanh(k*d) - w.^2;
        df = g*tanh(k*d) + g*k*d.*(sech(k*d).^2);
        k  = k - f./df;                   % converges to ~machine precision
    end
end
