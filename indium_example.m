%% Putting things together:

% Load data (optional afterwards)
% Calculate hook result

% Load hook result
% Recalculate change in heat loss estimate
% Get error in mass

% Calculate enthalpy
% Compare masses

%% Loads indium data from text files
% This block is optional after the first run — data can be loaded from
% a saved .mat file instead of re-reading all raw .txt files.

% clear
% dinfo = dir('*.txt');
% 
% % Pre-allocate arrays for heat flow (Q) and reference temperature (T)
% % Rows = time steps, Columns = 19 experimental runs
% indium_Q = zeros(10002,19);
% indium_T = zeros(10002,19);
% 
% for K = 1:19
%     fname = dinfo(K).name;
%     % Extract numeric portion of filename to determine run order
%     do = regexp(fname,'\d*','Match');
%     order(K) = str2double(do);
%     data = readmatrix(fname);
%     time_v = data(:,2);
%     tr_v = data(:,5);   % Reference temperature
%     ts_v = data(:,4);   % Sample temperature
%     q_v = data(:,3);    % Heat flow
%     stopper = length(time_v);
%     % Store only up to the actual length of this run (runs may differ in length)
%     indium_Q(1:stopper,K) = q_v;
%     indium_T(1:stopper,K) = tr_v;
%     indium_Ts(1:stopper,K) = ts_v;
% end

%% Load and concatenate data
% Load a single multi-segment DSC data file for indium calibration runs.
% The file spans -95 to 450°C at heating rates from 50 to 500 kK/s.
clear

% load indium_dat_set2.mat
% load indium_dat_set3.mat
% 
% % Data sets have different lengths (10002, 13008, 16674 rows).
% % Pad shorter arrays with zeros so they can be concatenated column-wise.
% indium_time(10003:16674,:) = 0;
% indium_Ts(10003:16674,:) = 0;
% indium_T(10003:16674,:) = 0;
% indium_Q(10003:16674,:) = 0;
% 
% in_Tr_set2(13009:16674,:) = 0;
% in_Ts_set2(13009:16674,:) = 0;
% in_t_set2(13009:16674,:) = 0;
% in_Q_set2(13009:16674,:) = 0;

%%
% Point to the specific indium calibration data file
dinfo = dir('cle_ind_baseline_06222026.txt');

for F = 1 : length(dinfo)  % Outer loop: iterate over matched files (typically just 1)
    fname = dinfo(F).name;  % Get filename
    % Read raw data, skipping the 5-line header
    data = readmatrix(fname, 'NumHeaderLines', 5);

    % Define column indices within the data matrix
    tsindex = 2; % Sample temperature (Ts), °C
    trindex = 3; % Reference temperature (Tr), °C
    vindex  = 4; % Heat flow, mW

    segment_number = 1;  % Tracks which DSC heating segment (ramp) we are in
    outputRow = 0;       % Row index within the current segment

    for inputRow = 1:size(data,1)  % Loop over every row in the raw data file

        outputRow = outputRow + 1;  % Advance output row counter

        % Segment boundary detection:
        % Rows with NaN in column 1 are separator lines — skip them
        if data(inputRow,1) == 'NaN'
            outputRow = outputRow - 1;  % Undo increment for this skipped row
            fprintf("test")
            continue
        % Rows with 0 in column 1 mark the start of a new segment
        elseif data(inputRow,1) == 0
            segment_number = segment_number + 1;  % Advance to next segment column
            outputRow = 1;                         % Reset row index for new segment
        end

        % Store current row into segment-indexed arrays.
        % Each column corresponds to one DSC heating ramp.
        t_all(outputRow,segment_number)  = data(inputRow,2);        % Time, s
        Ts_all(outputRow,segment_number) = data(inputRow,tsindex);  % Sample temp, °C
        Tr_all(outputRow,segment_number) = data(inputRow,trindex);  % Reference temp, °C (note: currently reads tsindex, not trindex)
        Q_all(outputRow,segment_number)  = data(inputRow,vindex);   % Heat flow, mW

    end

    % Compute the heating rate for each segment using two adjacent points
    % near row 25, where the ramp is expected to be well-established
    N = width(t_all);  % Total number of segments
    for k = 1:N
        HRs(k) = (Tr_all(25,k) - Tr_all(24,k)) / (t_all(25,k) - t_all(24,k));  % K/s
    end

    % Find the last valid (non-NaN) row in segment 2, used to trim all segments
    % to a consistent length for downstream concatenation
    [rowEnd] = find(isnan(t_all(:,2)), 1) - 1;

    % Extract segment 2 as the baseline heating ramp arrays
    time_heat        = t_all(1:rowEnd, 2);
    temp_sample_heat = Ts_all(1:rowEnd, 2);
    temp_ref_heat    = Tr_all(1:rowEnd, 2);
    heatflow_heat    = Q_all(1:rowEnd, 2);

    % Append segments 3 through end as additional columns (one per ramp)
    for ii = 3:segment_number
        time_heat        = cat(2, time_heat,        t_all(1:rowEnd, ii));
        temp_sample_heat = cat(2, temp_sample_heat, Ts_all(1:rowEnd, ii));
        temp_ref_heat    = cat(2, temp_ref_heat,    Tr_all(1:rowEnd, ii));
        heatflow_heat    = cat(2, heatflow_heat,    Q_all(1:rowEnd, ii));
    end


    %% Calculate peak areas
    % Find end of valid data in segment 2 (first NaN row)
    cend_1 = find(isnan(Q_all(:,2)), 1);

    for K = 2:N  % Loop over all segments
        % Find the last valid row for this segment
        cend = find(isnan(Q_all(:,K)), 1);
        if isempty(cend)
            cend = 35810;  % Fallback length if no NaN terminator is found
        end
        lengths_ind(K - 1) = cend - 1;  % Store usable length for each segment
        % Extract reference temperature and heat flow for this segment
        testT = Tr_all(1:cend-1, K)';  % Row vector of temperatures
        testQ = Q_all(1:cend-1, K)';   % Row vector of heat flows

        % Automatically detect the left and right baseline bounds of the
        % melting peak using the find_lefts_rights_mc function
        [T_ll, T_lr, T_rl, T_rr, ~, ~] = find_lefts_rights_mc(testT, testQ);

        % Integrate the melting peak area using a spline baseline correction.
        % Returns enthalpy (hcrys), corrected temperature (Tout), and heat flow (Qout).
        [hcrys, Tout, Qout] = spline_integral_orig(testT, testQ, T_ll, T_lr, T_rl, T_rr);

        % Store mean and std of enthalpy across the detected peak bounds
        hcrys_mean(K) = mean(hcrys);
        hcrys_std(K)  = std(hcrys);

        % Store mean onset and end temperatures for quality checking
        T_lr_means(K) = mean(T_lr);  % Mean left-right bound (onset side)
        T_rl_means(K) = mean(T_rl);  % Mean right-left bound (end side)

        % Reference values for indium:
        %   Enthalpy of fusion: 3283 J/mol (28.6 J/g)
        %   Cp at 300K: 26.98 J/K/mol (0.23 J/K/g)
        %   Source: https://akjournals.com/view/journals/10973/13/3/article-p419.xml

        %% Calculate hook heights and store values
        % The "hook" is a pre-peak deflection in the DSC signal caused by the
        % thermal lag between sample and sensor. Its height is used as an
        % independent mass estimator via the heat capacity (Cp).

        for P = 1:N
            if isnan(T_rl_means(:, 2))
            else
            testT = Tr_all(1:lengths_ind(:, 1),P);
            testQ = Q_all(1:lengths_ind(:, 1),P);
            excl = find(testT < T_lr_means(:, 2) | testT > T_rl_means(:, 2));
            try 
            [h_start, h_end, detail] = hook_MC_QAD(flipud(testT(excl)),flipud(testQ(excl)),0,1000);
            % h_start has form [value, noise]
            % h_end has form [mean value, std]
            % detail has form [mode, LOOCV SP, x_intersect (mean, std)]
            hookheights(P) = abs(h_end(1) - h_start(1));
            hstarts(P) = h_start(1);
            hends(P) = h_end(1);
            hooknoise(P) = h_start(2);
            hookstd(P) = h_end(2);
            x_inter_mean(P) = detail{3}(1);
            x_inter_std(P) = detail{3}(2);
            catch
                hookheights(P) = NaN;
                hooknoise(P) = NaN;
                hookstd(P) = NaN;
                x_inter_mean(P) = NaN;
                x_inter_std(P) = NaN;
            end
            end
        end
        

        %% Heat loss estimation (disabled — unreliable for large samples)
        % The correction below attempts to account for radiative/conductive heat
        % loss from the sample using the hook intercept temperature. It was found
        % to be unreliable for large samples and is left here for reference.
        % Indium Cp at 0°C: 0.23 J/g/K
        % for K = 1:19
        %     T_before = indium_Ts(1,K);
        %     Qloss(K) = hstarts(K)*(1-(x_inter_mean(K)./T_before));
        %     hook_mass(K) = (-hookheights(K)-Qloss(K))/HRs(K)/0.23/1e-6;
        % end

        %% Enthalpy peak using fixed peak bounds
        % As a consistency check, re-integrate each peak using fixed temperature
        % bounds (125–140°C onset side, 165–180°C end side) rather than the
        % automatically detected bounds above
        for K = 1:19
            [hcrys, ~, ~] = spline_integral_orig(testT, testQ, 125, 140, 165, 180);
            hcrys_fix(K) = hcrys(1);  % Store only the first (primary) integral value
        end
    end

    %% Calculate all the masses
    % Estimate sample mass from enthalpy method: divide integrated heat (mJ)
    % by heating rate (K/s) and indium's specific enthalpy of fusion (28.6 J/g),
    % converting to nanograms (factor of 1e-9)
    mass_hcrys = abs(hcrys_mean * 1e-3 ./ HRs ./ 28.6 / 1e-9);

    % Estimate sample mass from hook method: divide hook height (mW offset)
    % by heating rate and indium's Cp (0.23 J/g/K), converting to nanograms
    mass_hook = abs(hookheights * 1e-3 ./ HRs ./ 0.23 / 1e-9);

    %% Compare enthalpies (runs 1–9): mass estimates vs. heating rate
    figure()
    % Enthalpy-derived mass (moving peak bounds) — filled blue dots
    semilogx(HRs(1:9), abs(hcrys_mean(1:9)) * 1e-3 ./ Q_all(1:9) ./ 28.6 / 1e-9, ...
        '.b', 'markersize', 36, 'linewidth', 2)
    hold on
    % Hook-derived mass — blue crosses
    semilogx(HRs(1:9), abs(hookheights(1:9)) * 1e-3 ./ Q_all(1:9) ./ 0.23 / 1e-9, ...
        'xb', 'markersize', 12, 'linewidth', 2)
    xlim([4e-1 2e4])
    ylim([1e3 6e3])
    ylabel("Estimated mass (ng)")
    xlabel("Heating rate (K/s)")
    set(gca, 'fontsize', 24)
    legend("Enthalpy (moving peak)", "Enthalpy (fixed peak)", "Hook", "Hook with heat loss", ...
        'location', 'best')

    % %% Compare enthalpy vs. hook mass estimates as a function of sample mass
    % % Restrict to runs with physically meaningful heating rates (0.1–2000 K/s)
    excl2 = find(Q_all < 2e3 & Q_all > 0.1);
    dummy = linspace(5e-3, max(mass_hook(1:9)) + 500);  % 1:1 reference line range

    figure()
    loglog(mass_hcrys(1:9), mass_hook(1:9), 'or', 'markersize', 12)
    hold on
    loglog(dummy, dummy)  % 1:1 line — perfect agreement between both mass estimators
    % Note: agreement is sensitive to accuracy of the heating rate estimate

    %% Ratio of hook mass to enthalpy mass as a function of heating rate
    % A ratio of 1.0 indicates the two methods agree; shaded bands show ±10% and ±20%
    valid = 2:N;  % indices where hookheights and hcrys_mean are actually populated
    subs = valid(Q_all(valid) > 0.5);  % filter by heating rate within valid range
    figure()
    semilogx(NaN, NaN)  % Dummy call to initialize axes before patch() calls

    % Shaded band: ±10% agreement (darker blue)
    patch([5e-1 2e4 2e4 5e-1], [0.9 0.9 1.1 1.1], [111, 175, 237]/256, 'edgecolor', 'none')
    hold on
    % Shaded band: 10–20% over-estimate (lighter blue, upper)
    patch([5e-1 2e4 2e4 5e-1], [1.1 1.1 1.2 1.2], [183, 216, 247]/256, 'edgecolor', 'none')
    % Shaded band: 10–20% under-estimate (lighter blue, lower)
    patch([5e-1 2e4 2e4 5e-1], [0.8 0.8 0.9 0.9], [183, 216, 247]/256, 'edgecolor', 'none')

    % Plot ratio: hook mass / enthalpy mass for each valid segment
    semilogx(HRs(subs), ...
        (abs(hookheights(subs)) * 1e-3 ./ Q_all(subs) ./ 0.23 / 1e-9) ./ ...
        (abs(hcrys_mean(subs)) * 1e-3 ./ Q_all(subs) ./ 28.6 / 1e-9), ...
        '.b', 'markersize', 36, 'linewidth', 2)

    yline(1, 'k:', 'linewidth', 2)  % Reference line at ratio = 1 (perfect agreement)
    ylabel("Ratio between mass estimates")
    xlabel("Heating rate (K/s)")
    set(gca, 'fontsize', 24)
    xlim([5e-1 2e4])
    ylim([0.75 1.25])

end