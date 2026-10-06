function [hcrys, Tout, Qout] = spline_integral_orig(T_data, Q_data, T_ll, T_lr, T_rl, T_rr)
% spline_integral_orig  Integrate a DSC melting peak using a spline baseline correction.
%
% For each of the 1000 Monte Carlo boundary samples, this function:
%   1. Fits a linear baseline to the left and right baseline regions
%   2. Interpolates a spline through those two linear fits
%   3. Subtracts the spline baseline from the raw heat flow
%   4. Integrates the corrected peak area via the trapezoidal rule
%
% Inputs:
%   T_data  - Temperature vector (°C)          
%   Q_data  - Heat flow vector (mW)            
%   T_ll    - 1000x1 vector, left outer baseline bound temperatures (°C)
%   T_lr    - 1000x1 vector, left inner (peak onset) bound temperatures (°C)
%   T_rl    - 1000x1 vector, right inner (peak end) bound temperatures (°C)
%   T_rr    - 1000x1 vector, right outer baseline bound temperatures (°C)
%   figflag - Plot flag (1 = show diagnostic figure, 0 = suppress)
%             [NOTE: plotting block is currently commented out]
%
% Outputs:
%   hcrys - 1x1000 vector of integrated peak enthalpies (mW·K per sample),
%           one value per Monte Carlo boundary draw
%   Tout  - Matrix of baseline-corrected temperature vectors (one column per MC draw)
%   Qout  - Matrix of baseline-corrected heat flow vectors (one column per MC draw)
%

    % --- Initialise output array ---
    % hcrys stores one integrated enthalpy value per Monte Carlo boundary draw
    % Nov. 2024 version — accepts vectors for each bound point, one column of T/Q at a time

    hcrys = zeros(length(T_ll), 1);  % 1 x number_of_segments


    T_current = T_data(:);
    Q_current = Q_data(:);

    % Pre-allocated output matrices 
    Tout = zeros(length(T_data), length(T_data(1,:)));
    Qout = zeros(length(T_data), length(T_data(1,:)));

    % --- Monte Carlo integration loop ---
    % Iterate over all 1000 boundary sample sets
    for c = 2:length(T_ll)

        % Find indices of data points falling within each boundary window
        left      = find(T_current > T_ll(c) & T_current < T_lr(c));  % Left baseline region
        right     = find(T_current > T_rl(c) & T_current < T_rr(c));  % Right baseline region
        corrected = find(T_current > T_ll(c) & T_current < T_rr(c));  % Full peak + baseline window
        % Extract temperature and heat flow for each region
        T_left      = T_current(left);
        T_right     = T_current(right);
        Q_left      = Q_current(left);
        Q_right     = Q_current(right);
        T_corrected = T_current(corrected);  % Temperature over the full integration window
        T_corrected = T_corrected(:);
        Q_temp      = Q_current(corrected);  % Raw heat flow over the full integration window
        Q_temp = Q_temp(:);
        %force into column vectors
        T_right = T_right(:);
        Q_right = Q_right(:);
        
        T_left  = T_left(:);
        Q_left  = Q_left(:);

        % --- Fit linear baselines to the left and right anchor regions ---
        % Using least-squares: [ones, T] \ Q  gives [intercept; slope]

        fit_right = [ones(length(T_right), 1) T_right] \ (Q_right - 1);  % Right baseline: [b; m]
        fit_left  = [ones(length(T_left),  1) T_left]  \ (Q_left - 1);   % Left baseline:  [b; m]


        % Evaluate fitted lines at their respective temperature points
        T_fit = [T_right; T_left];  % Combined temperature knots for spline
        Q_fit = [fit_right(1) + fit_right(2)*T_right; ...
                 fit_left(1)  + fit_left(2) *T_left];  % Corresponding baseline values

        % Remove any duplicate temperature values before spline fitting
        % (duplicate knots would cause spline() to error)
        [T_fit, uidT, ~] = unique(T_fit);
        Q_fit = Q_fit(uidT);

        % --- Build and evaluate a spline baseline across the full peak window ---
        % The spline smoothly connects the left and right linear fits,
        % providing a continuous baseline under the melting peak
        try
            Q_spline = spline(T_fit, Q_fit, T_corrected);
            % if ~isvector(Q_spline)
            %     fprintf('Iteration c=%d: Q_spline is %dx%d\n', c, size(Q_spline,1), size(Q_spline,2))
            %     fprintf('  T_fit: %dx%d, Q_fit: %dx%d, T_corrected: %dx%d\n', ...
            %     size(T_fit,1), size(T_fit,2), ...
            %     size(Q_fit,1), size(Q_fit,2), ...
            %     size(T_corrected,1), size(T_corrected,2))
            %  end
            Q_spline = Q_spline(:);
        catch
            % If spline fails, print sizes to help diagnose the mismatch
            disp(size(T_fit))
            disp(size(Q_fit))
            disp(size(T_corrected))
        end

        % --- Baseline correction and integration ---
        % Subtract the spline baseline from the raw heat flow
        
        Q_corr = abs(Q_temp - Q_spline);
        
        % Integrate the corrected peak using the trapezoidal rule.
        % Result is in mW·K, proportional to the enthalpy of fusion.
        % --- Baseline correction and integration ---
        
        % disp('hcrys'); disp(size(hcrys))
        %     disp('T_corrected'); disp(size(T_corrected))
        %     disp('Q_corr'); disp(size(Q_corr))
        % Only integrate if there are enough points to form a valid integral
        if length(T_corrected) > 1 && length(Q_corr) > 1
            hcrys(c) = trapz(T_corrected, Q_corr);
            
        else
            hcrys(c) = NaN;  % Flag this draw as invalid rather than crashing
        end


        % Store the corrected temperature and heat flow for this MC draw.
        % Columns correspond to MC sample index c.
        Tout(1:length(T_corrected), c) = T_corrected;
        Qout(1:length(Q_corr),      c) = Q_corr;

        % --- Diagnostic plot (disabled) ---
        % Uncomment figflag block to visualise baseline correction for each draw:
        %   p1 = T_fit
        %   bounds = [T_ll T_lr T_rl T_rr]
        %   figure()
        %   plot(T_data, Q_data)
        %
        % if figflag
        %     figure();
        %     hold on
        %     plot(T_data, Q_data)
        %     plot(T_corrected, Q_temp,   'b',  'LineWidth', 0.5)  % Raw data in window
        %     plot(T_fit,       Q_fit,    'b.', 'LineWidth', 1)    % Baseline knots
        %     plot(T_corrected, Q_spline, 'b--','LineWidth', 0.5)  % Spline baseline
        %     plot(T_corrected, Q_temp,   'c',  'LineWidth', 0.5)  % Raw (duplicate, cyan)
        %     plot(T_fit,       Q_fit,    'c.', 'LineWidth', 1)    % Knots (duplicate, cyan)
        %     plot(T_corrected, Q_spline, 'c--','LineWidth', 0.5)  % Spline (duplicate, cyan)
        %     plot(T_corrected, Q_corr)                            % Corrected peak
        %     xline(T_lr)   % Peak onset boundary
        %     xline(T_rl)   % Peak end boundary
        %     yline(0)      % Zero heat flow reference
        %     hold off
        % end

    end  % End Monte Carlo loop

end  % End file-reading loop