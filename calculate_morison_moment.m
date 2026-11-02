function M = calculate_morison_moment(comp, t, eta, geom, hydro, z_weld, d, opt)
%CALCULATE_MORISON_MOMENT  Wave bending moment at one or many welds (Airy + Morison).
%
%   M = CALCULATE_MORISON_MOMENT(comp, t, eta, geom, hydro, z_weld, d, opt)
%
%   z_weld may be a SCALAR or a VECTOR of Nw weld elevations [mMSL]. The wave
%   kinematics and the distributed Morison force f(z,t) [N/m] are built ONCE on a
%   shared vertical grid spanning from the lowest weld up to just above the crest;
%   the bending moment is then obtained for EACH weld by integrating
%   f(z,t)*(z - z_weld(i)) from that weld up to the instantaneous surface eta(t).
%   The expensive O(Nz*Nfreq*Nt) kinematics matmul is therefore paid ONCE per sea
%   state regardless of weld count -- a 37-weld run costs essentially the same as
%   one weld plus 37 cheap masked-trapz integrations. (This is the key efficiency
%   requirement: kinematics are NOT recomputed per weld.)
%
%   KINEMATICS IN THE SPLASH ZONE (opt.kin) -- DNV-RP-C205 Sec.3.3.3
%     'extrapolated' (DEFAULT) : linear velocity profile held constant at its
%                                still-water-level value for 0 < z < eta (cosh
%                                argument clamped, zc = min(z,0)); force integrated
%                                to eta. DNV's recommendation for a LINEAR input
%                                surface elevation (a random-phase sum of cosines).
%     'still_water'            : integrate only to z = 0 (ignore the splash zone).
%     'wheeler'                : Wheeler (1970) stretching. DNV-RP-C205 [3.3.3.4]
%                                advises AGAINST it for a linear input elevation
%                                (it underpredicts kinematics around/below z = 0);
%                                a warning is raised if selected.
%
%   INPUTS
%     comp    struct from generate_jonswap_timeseries (.w .k .a .phi, each Nfreq x 1)
%     t       time vector              [s]   (1 x Nt)
%     eta     surface elevation        [m]   (1 x Nt)
%     geom    .ODfun(z) steel OD [m], .mgfun(z) marine-growth thickness [m]
%     hydro   .rho [kg/m^3], .Cd [-], .Cm [-], .current(z) [m/s]
%     z_weld  weld elevation(s)        [mMSL]  (scalar OR Nw x 1; e.g. -d at mudline)
%     d       mean water depth (MSL -> mudline) [m]
%     opt     .Nz nodes, .kin model (default 'extrapolated'), .z_top [mMSL] (opt),
%             .dt [s], .applyDAF [logical], .fn [Hz], .zeta [-]
%
%   OUTPUT
%     M       bending moment at the weld(s)  [N*m]  (Nw x Nt) -- one row per weld
%
%   Conventions: z is elevation in mMSL (0 = MSL, -d = mudline, +ve up); pile at
%   x = 0; in-line (wave-direction) moment only (omni-directional, no misalignment).
%   Requires MATLAB R2016b+ (implicit expansion).
% -------------------------------------------------------------------------

    z_weld = z_weld(:);                          % Nw x 1
    Nw     = numel(z_weld);
    Nt     = numel(t);

    w   = comp.w;  k = comp.k;  a = comp.a;  phi = comp.phi;   % Nfreq x 1
    rho = hydro.rho;  Cd = hydro.Cd;  Cm = hydro.Cm;
    Nz  = opt.Nz;
    if isfield(opt,'kin') && ~isempty(opt.kin), kin = lower(opt.kin); else, kin = 'extrapolated'; end

    % --- phase / trig: largest arrays, built ONCE (Nfreq x Nt) ------------
    Theta = w*t + phi;
    C     = cos(Theta);
    Sn    = sin(Theta);
    clear Theta
    AW     = a.*w;                               % Nfreq x 1
    AW2    = a.*w.^2;                            % Nfreq x 1
    sinhkd = sinh(k*d);                          % Nfreq x 1

    M = zeros(Nw, Nt);

    switch kin
    % =====================================================================
    case {'extrapolated','still_water'}
        % shared PHYSICAL grid: lowest weld -> just above the highest crest
        if isfield(opt,'z_top') && ~isempty(opt.z_top), z_top = opt.z_top; else, z_top = max(eta) + 2; end
        z_floor = min(z_weld);
        if z_top <= z_floor
            error('calculate_morison_moment:ztop','z_top (%.2f) must exceed the lowest weld (%.2f).', z_top, z_floor);
        end
        z   = linspace(z_floor, z_top, Nz).';     % Nz x 1  (spans all welds)
        zc  = min(z, 0);                          % clamp profile at MSL
        Phi = cosh((zc + d) * k.') ./ (sinhkd.'); % Nz x Nfreq   (built once)

        % kinematics + force on the full column, ONCE (Nz x Nt) --------------
        U  =  (Phi .* AW.')  * C;                 % horizontal velocity
        AC = -(Phi .* AW2.') * Sn;                % horizontal acceleration
        clear C Sn
        Dh   = geom.ODfun(z) + 2*geom.mgfun(z);   % Nz x 1  hydrodynamic diameter
        Uc   = hydro.current(z);                  % Nz x 1
        Utot = U + Uc;
        f = 0.5*rho*Cd .* Dh .* abs(Utot).*Utot ...   % drag
            + rho*Cm .* (pi*Dh.^2/4) .* AC;           % inertia   (Nz x Nt)
        clear U AC Utot

        if strcmp(kin,'still_water')
            surf = (z <= 0);                      % Nz x 1  integrate to MSL
        else
            surf = (z <= eta);                    % Nz x Nt integrate to surface
            if max(eta) > z_top
                warning('calculate_morison_moment:ztop', ...
                    'max(eta)=%.2f m exceeds z_top=%.2f m; raise opt.z_top.', max(eta), z_top);
            end
        end

        % per-weld moment: cheap masked trapz (kinematics already done) ------
        for i = 1:Nw
            wet   = (z >= z_weld(i)) & surf;      % Nz x 1 or Nz x Nt
            lever = (z - z_weld(i));              % Nz x 1
            M(i,:) = trapz(z, f .* lever .* wet, 1);
        end

    % =====================================================================
    case 'wheeler'
        warning('calculate_morison_moment:wheeler', ...
           ['Wheeler stretching with a LINEAR input elevation is advised against by ', ...
            'DNV-RP-C205 [3.3.3.4]: it underpredicts kinematics near and below z = 0. ', ...
            'Use opt.kin = ''extrapolated'' for a random-phase (linear) sea.']);

        zeta = linspace(-d, 0, Nz).';             % computational grid
        Phi  = cosh((zeta + d) * k.') ./ (sinhkd.');  % Nz x Nfreq (once)
        U  =  (Phi .* AW.')  * C;                 % Nz x Nt (once)
        AC = -(Phi .* AW2.') * Sn;
        clear C Sn
        Z    = zeta * ((d + eta)/d) + eta;        % Nz x Nt physical elevations
        Dh   = geom.ODfun(Z) + 2*geom.mgfun(Z);   % Nz x Nt
        Uc   = hydro.current(Z);                  % Nz x Nt
        Utot = U + Uc;
        f = 0.5*rho*Cd .* Dh .* abs(Utot).*Utot ...
            + rho*Cm .* (pi*Dh.^2/4) .* AC;       % Nz x Nt
        clear U AC Utot Dh Uc
        Jac = ((d + eta)/d);                      % 1 x Nt  (uniform-zeta Jacobian)

        for i = 1:Nw
            wet   = (Z >= z_weld(i));             % Nz x Nt
            lever = Z - z_weld(i);                % Nz x Nt
            M(i,:) = trapz(zeta, f .* lever .* wet, 1) .* Jac;
        end

    otherwise
        error('calculate_morison_moment:kin', ...
            'Unknown opt.kin = "%s" (use extrapolated | still_water | wheeler).', kin);
    end

    % --- optional SDOF dynamic amplification (per weld row) ---------------
    if isfield(opt,'applyDAF') && opt.applyDAF
        for i = 1:Nw
            M(i,:) = applySDOFdaf(M(i,:), opt.dt, opt.fn, opt.zeta);
        end
    end
end

% =========================================================================
function Mdyn = applySDOFdaf(M, dt, fn, zeta)
%APPLYSDOFDAF  Frequency-domain SDOF dynamic amplification of a moment series.
%   Mdyn(t) = IFFT{ H(f) * FFT[M(t)] },  H = 1/((1-r^2) + i 2 zeta r),  r = f/fn.
    M   = M(:).';
    N   = numel(M);
    f   = (0:N-1)/(N*dt);
    f(f > 1/(2*dt)) = f(f > 1/(2*dt)) - 1/dt;
    r   = f/fn;
    H   = 1 ./ ((1 - r.^2) + 1i*2*zeta.*r);
    Mdyn = real(ifft(fft(M).*H));
end
