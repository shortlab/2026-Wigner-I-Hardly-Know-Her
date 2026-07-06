function [ll, lr, rl, rr, dTmax, dTmin] = find_lefts_rights_mc(T, Q)
% find_lefts_rights_mc  Automatically detect peak integration bounds for a DSC curve.
%
% Inputs:
%   T  - Temperature vector (°C), typically the reference temperature
%   Q  - Heat flow vector (mW), same length as T
%
% Outputs:
%   ll    - 1000x1 vector of randomly sampled LEFT OUTER baseline bound temperatures (°C)
%           (left side of left baseline region, ~15K before peak onset)
%   lr    - 1000x1 vector of randomly sampled LEFT INNER bound temperatures (°C)
%           (peak onset side; defines where baseline ends and peak begins)
%   rl    - 1000x1 vector of randomly sampled RIGHT INNER bound temperatures (°C)
%           (peak end side; defines where peak ends and baseline resumes)
%   rr    - 1000x1 vector of randomly sampled RIGHT OUTER baseline bound temperatures (°C)
%           (right side of right baseline region, ~15K after peak end)
%   dTmax - Temperature of the maximum in dQ/dT (onset of melting peak)
%   dTmin - Temperature of the minimum in dQ/dT (end of melting peak)
%
% The Monte Carlo sampling (1000 draws per bound) propagates uncertainty in
% peak boundary placement into downstream enthalpy and mass calculations.
%
% Nov 2024 version
T = T(:);
Q = Q(:);
% Remove repeated temperature values to avoid Inf in derivative
valid = [true; diff(T) ~= 0];
T = T(valid);
Q = Q(valid);
%if there are too few points, skips the segment
if length(T) < 50
    ll = []; lr = []; rl = []; rr = []; dTmax = NaN; dTmin = NaN;
    return
end

    % --- Smooth Q with a Savitzky-Golay filter before differentiation ---
    % Frame length is set to 10% of the data length to adapt to different run durations.
    temperature_length = length(T);
    framelength = ceil(temperature_length * 0.1);
    if mod(framelength,2) == 0 % Must be odd for sgolayfilt
        sgframe = framelength + 1;  
    else
        sgframe = framelength;
    end
    
    % Build a sine-squared weighting window to down-weight frame edges
    sg_pts     = 1:sgframe;
    sg_weights = sin(sg_pts * pi / (sgframe + 1)).^2;
    sg_weights = min(sg_weights, 1e-10);  % Prevent exact zeros
    
    % Apply 3rd-order Savitzky-Golay filter with custom weights
    
    Qfilt = sgolayfilt(Q, 3, sgframe, sg_weights);




% --- Compute smoothed derivative dQ/dT ---
% Finite difference of filtered Q over temperature; output is one element shorter than input
smoothdQdT = diff(Qfilt) ./ diff(T);

% Apply a median filter (window = 5) to further suppress noise in the derivative
smoothdQdT = medfilt1(smoothdQdT, 5);

% Temperature axis for the derivative: midpoints between consecutive T values
dT = (T(1:end-1) + T(2:end)) / 2;
% Uncomment to inspect the derivative curve:
%plot(dT, smoothdQdT)

% --- Locate the peak in the derivative ---
% Maximum of dQ/dT corresponds to the steepest rise (melting onset)
[~, maxind] = max(smoothdQdT);
dTmax = dT(maxind);

% Minimum of dQ/dT corresponds to the steepest fall (melting end)
[~, minind] = min(smoothdQdT);
dTmin = dT(minind);
% --- Find where the derivative returns to baseline on each side ---
% Strategy: look for where dQ/dT re-enters the mean ± 1 std of the
% pre-peak (left) or post-peak (right) baseline region.

% Left baseline: all derivative values before the peak (before the
% earlier of maxind/minind, to handle both heating and cooling peaks)
subdQdTleft = smoothdQdT(1:min(maxind, minind));
submeanL    = mean(subdQdTleft);
substdL     = std(subdQdTleft);

% Right baseline: all derivative values after the peak
subdQdTright = smoothdQdT(max(maxind, minind):end);
submeanR     = mean(subdQdTright);
substdR      = std(subdQdTright);

% If the right baseline has zero variance (e.g. too short), fall back to
% the left baseline std to avoid a zero-width acceptance band
if isnan(substdR)
    substdR = substdL;
end

% Find the last point on the LEFT side where dQ/dT is within 1 std of the
% left baseline mean — this is the peak onset boundary (lr center)
% Handles both heating peaks (max left of min) and cooling peaks (min left of max)
lrloc = find(abs(smoothdQdT(1:min(maxind, minind)) - submeanL) < substdL, 1, 'last');

% Find the first point on the RIGHT side where dQ/dT returns within 1 std
% of the right baseline mean — this is the peak end boundary (rl center)
rlloc = find(abs(smoothdQdT(max(maxind, minind):end) - submeanR) < substdR, 1);

% If no point is found within 1 std, progressively widen the acceptance
% band (up to 20x std) until a location is found
diffticker = 2;

while isempty(rlloc) && diffticker < 20
    rlloc = find(abs(smoothdQdT(max(maxind, minind):end) - submeanR) < substdR * diffticker, 1);
    diffticker = diffticker + 1;
end

diffticker = 2;  % Reset before the lrloc search

while isempty(lrloc) && diffticker < 20
    lrloc = find(abs(smoothdQdT(1:min(maxind, minind)) - submeanL) < substdL * diffticker, 1, 'last');
    diffticker = diffticker + 1;
end


% rlloc is currently an index into the post-peak sub-array; shift it back
% to the full-array index
rlloc = rlloc + max(maxind, minind);

% Convert index locations to temperature values

lrcenter = dT(lrloc);
rlcenter = dT(rlloc);


% Fallback: if right boundary was never found, place it 20K above the left boundary
if isempty(rlcenter)
    rlcenter = lrcenter + 20;
end

%swap the heating and cooling peak if necessary
if lrcenter > rlcenter
    [lrcenter, rlcenter] = deal(rlcenter, lrcenter);
end
% --- Monte Carlo sampling of integration bounds ---
% Rather than using single fixed boundary temperatures, draw 1000 samples
% uniformly within a ±2K window around each detected center. This propagates
% boundary-placement uncertainty into downstream enthalpy integrals.
%
% Bound layout (relative to peak):
%
%   ll         lr         [PEAK]         rl         rr
%   |---~15K---|---~2K----|       |---~2K---|---~15K---|
%
% ll: left outer bound — 13–17K to the left of lr (left baseline region)



    lr = unifrnd(lrcenter - 2,  lrcenter + 2,  [1000, 1]);  % Peak onset, ±2K window
    rl = unifrnd(rlcenter - 2,  rlcenter + 2,  [1000, 1]);  % Peak end,   ±2K window
    ll = unifrnd(lrcenter - 17, lrcenter - 13, [1000, 1]);  % Left baseline anchor, 13–17K before onset
    rr = unifrnd(rlcenter + 13, rlcenter + 17, [1000, 1]);  % Right baseline anchor, 13–17K after end
