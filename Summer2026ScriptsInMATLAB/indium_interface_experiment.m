clear all
dinfo = dir('*.txt');

% Cell arrays to store data from each file
all_temp = {};
all_heatflow = {};

for K = 1 : length(dinfo) %looks through all the files (1)
    fname = dinfo(K).name;%smush into a matrix
    data = readmatrix(fname, 'NumHeaderLines', 5);
    tsindex = 2; % Ts
    trindex = 3; % Tr
    vindex  = 4; % mW
    segment_number=1;
    outputRow=0;

    % Clear segment arrays for each file
    clear time_all temp_sample_all temp_ref_all heatflow_all

for inputRow=1:size(data,1)%for all the rows, do this
    outputRow=outputRow+1;%incrementing
% Split into each segment: looking for the NaN endings
if data(inputRow,1) == 'NaN'
        outputRow=outputRow-1;
        fprintf("test")
continue
elseif data(inputRow,1) == 0
        segment_number=segment_number+1;
        outputRow=1;
end
% Allocate data into variables
    time_all(outputRow,segment_number)=data(inputRow,2);           % s
    temp_sample_all(outputRow,segment_number)=data(inputRow,tsindex);    % C
    temp_ref_all(outputRow,segment_number)=data(inputRow,tsindex);       % C
    heatflow_all(outputRow,segment_number)=data(inputRow,vindex);       % mW
end
[rowEnd]=find(isnan(time_all(:,2)),1)-1;
time_heat=time_all(1:rowEnd,2);
temp_sample_heat=temp_sample_all(1:rowEnd,2);
temp_ref_heat=temp_ref_all(1:rowEnd,2);
heatflow_heat=heatflow_all(1:rowEnd,2);
for ii=3:segment_number
time_heat=cat(2,time_heat,time_all(1:rowEnd,ii));
temp_sample_heat=cat(2,temp_sample_heat,temp_sample_all(1:rowEnd,ii));
temp_ref_heat=cat(2,temp_ref_heat,temp_ref_all(1:rowEnd,ii));
heatflow_heat=cat(2,heatflow_heat,heatflow_all(1:rowEnd,ii));
end

    % Store this file's data
    all_temp{K} = temp_sample_heat;
    all_heatflow{K} = heatflow_heat;
end

%% Separate figure per file — subtract indium baseline only for "ind" files
%% (e.g. copper oxide files) so their data is shown net of the indium signal
%% All figures share the same y-axis scale

% Find the baseline file index by name
baseline_idx = find(contains({dinfo.name}, 'cle_ind_baseline'));
if isempty(baseline_idx)
    error('Could not find a file named cle_ind_baseline in dinfo.');
elseif numel(baseline_idx) > 1
    error('Multiple files match cle_ind_baseline — name collision.');
end

temp_base = all_temp{baseline_idx};
flow_base = all_heatflow{baseline_idx};

%% First pass: compute the y-data for every file/segment and track global min/max
plot_data = cell(length(dinfo), 1);  % store {seg}.x, {seg}.y per file
global_ymin = Inf;
global_ymax = -Inf;

for fidx = 1:length(dinfo)
    temp1 = all_temp{fidx};
    flow1 = all_heatflow{fidx};
    numSegs = size(temp1,2);

    % Files with 'ind' in the name (e.g. copper oxide runs on the indium
    % chip) get the indium baseline subtracted; the baseline file itself
    % is shown raw.
    is_ind_file = contains(dinfo(fidx).name, 'ind') && fidx ~= baseline_idx;

    seg_x = cell(numSegs,1);
    seg_y = cell(numSegs,1);

    for seg = 1:numSegs
        if is_ind_file
            numSegsUse = min(numSegs, size(temp_base,2));
            if seg > numSegsUse
                continue
            end
            flow_base_interp = interp1(temp_base(:,seg), flow_base(:,seg), temp1(:,seg), 'linear', NaN);
            y = flow1(:,seg) - flow_base_interp;
        else
            y = flow1(:,seg);
        end
        seg_x{seg} = temp1(:,seg);
        seg_y{seg} = y;

        global_ymin = min(global_ymin, min(y, [], 'omitnan'));
        global_ymax = max(global_ymax, max(y, [], 'omitnan'));
    end

    plot_data{fidx} = struct('x', {seg_x}, 'y', {seg_y}, 'is_ind_file', is_ind_file);
end

% Add a little padding so curves don't touch the plot edges
yrange = global_ymax - global_ymin;
ypad = 0.05 * yrange;
ylims = [global_ymin - ypad, global_ymax + ypad];

%% Second pass: plot each file with the shared y-axis scale
for fidx = 1:length(dinfo)
    data = plot_data{fidx};
    numSegs = numel(data.x);

    figure(fidx)
    clf   % clear the figure so old plots don't linger between runs
    hold on
    for seg = 1:numSegs
        if isempty(data.x{seg})
            continue
        end
        plot(data.x{seg}, data.y{seg}, '-', 'LineWidth', 0.5, 'DisplayName', num2str(seg))
    end
    hold off

    set(gca,'FontSize',15)
    xlabel('Temperature (\circC)','FontSize',20)
    if data.is_ind_file
        ylabel('\Delta Power (mW)','FontSize',20)
        title(sprintf('%s (indium baseline subtracted)', dinfo(fidx).name), 'Interpreter','none')
        yline(0, '--k', 'LineWidth', 0.5)
    else
        ylabel('Power (mW)','FontSize',20)
        title(dinfo(fidx).name, 'Interpreter','none')
    end
    xlim([20 325])
    ylim(ylims)
    legend('show')
end