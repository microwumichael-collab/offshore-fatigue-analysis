function scatter_data = build_scatter_jonswap(g, gamma_js, db_extremes)
%BUILD_SCATTER_JONSWAP  Build fatigue scatter from main_fatigue_analysis constants.
%
%   scatter_data = BUILD_SCATTER_JONSWAP(g, gamma_js, db_extremes)
%
%   INPUTS  (taken directly from main_fatigue_analysis.m)
%     g            gravity [m/s^2]                        (9.81)
%     gamma_js     JONSWAP peak-enhancement factor [-]    (3.3)
%     db_extremes  3x2 matrix [Hs, Tp] from DB Table 12-2
%                  e.g. [6.3,10.4; 11.0,12.9; 13.5,14.0]
%
%   OUTPUT
%     scatter_data  Nstate x 3 matrix  [Hs, Tp, prob]
%                   Drop-in replacement for the scatter_data block.

    % --- input guard ------------------------------------------------------
    assert(ismatrix(db_extremes) && size(db_extremes,2) == 2, ...
        'db_extremes must be an N x 2 matrix of [Hs, Tp] pairs.');

    % --- Hs bins: 0.5 m spacing, capped just above 1-yr return Hs --------
    Hs_bins = (0.5 : 0.5 : 7.0).';         % 14 bins [m]

    % =====================================================================
    % STEP 1  Calibrate JONSWAP peak steepness from DB extreme pairs
    %         sp = 2*pi*Hs / (g*Tp^2)   [DNV-RP-C205 Sec. 3.5.5]
    % =====================================================================
    Hs_ext   = db_extremes(:,1);            % column 1: Hs values
    Tp_ext   = db_extremes(:,2);            % column 2: Tp values
    sp_each  = 2*pi * Hs_ext ./ (g * Tp_ext.^2);
    sp       = mean(sp_each);

    % Tp for each bin from calibrated steepness
    Tp_bins  = sqrt(2*pi * Hs_bins / (g * sp));

    % =====================================================================
    % STEP 2  Integrate JONSWAP spectrum -> m0 -> occurrence weights
    %         prob ∝ 1/m0 = 16/Hs^2   [DNV-RP-C205 Sec. 3.5.3]
    % =====================================================================
    Nw_int = 4096;
    m0_v   = zeros(size(Hs_bins));

    for k = 1:numel(Hs_bins)
        Hs_k = Hs_bins(k);
        wp_k = 2*pi / Tp_bins(k);

        w               = linspace(0.05, 4.0*wp_k, Nw_int);
        dw              = w(2) - w(1);
        sigma           = 0.07 * ones(size(w));
        sigma(w > wp_k) = 0.09;
        Ag  = 1 - 0.287 * log(gamma_js);
        Spm = (5/16) * Hs_k^2 * wp_k^4 .* w.^(-5) .* ...
              exp(-1.25 * (w/wp_k).^(-4));
        S   = Ag * Spm .* ...
              gamma_js.^(exp(-0.5*((w - wp_k)./(sigma*wp_k)).^2));

        m0_v(k) = sum(S) * dw;
    end

    prob_bins = (1./m0_v) / sum(1./m0_v);  % normalised, sum = 1

    % =====================================================================
    % Output
    % =====================================================================
    scatter_data = [Hs_bins, Tp_bins, prob_bins];

    % console summary
    fprintf('\n  build_scatter_jonswap: sp=%.4f  Tp=%.3f*sqrt(Hs)\n', ...
            sp, sqrt(2*pi/(g*sp)));
    fprintf('  %-6s  %-7s  %-8s\n','Hs[m]','Tp[s]','prob[-]');
    fprintf('  %s\n', repmat('-',1,26));
    for k = 1:numel(Hs_bins)
        fprintf('  %5.1f   %6.2f   %7.5f\n', ...
                scatter_data(k,1), scatter_data(k,2), scatter_data(k,3));
    end
    fprintf('  sum(prob) = %.6f\n\n', sum(prob_bins));

end