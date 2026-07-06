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

%% Plot raw data for each heating segment
figure(K)
% splitting into heating vs cooling segments
flow_heat = heatflow_heat(:,3:2:end-2);
t_heat = temp_sample_heat(:,3:2:end-2);
flow_cool = heatflow_heat(:,2:2:end-2);
t_cool = temp_sample_heat(:,2:2:end-2);
for heat=1:segment_number-1
    plot(temp_sample_heat(:,heat),heatflow_heat(:,heat),'-','LineWidth',0.5,'DisplayName',num2str(heat))
    hold on
end
set(gca,'FontSize',15)
xlabel('Temperature (\circC)','FontSize',20)
ylabel('Power (mW)','FontSize',20)
title(fname,'Interpreter','none')
xlim([20 325])
hold on
% legend('Location','best')
end

%% Plot file2 - file1 (subtracted) if at least 2 files exist
    temp1 = all_temp{2};
    temp2 = all_temp{3};
    flow1 = all_heatflow{2};
    flow2 = all_heatflow{3};

   
    numSegs = min(size(temp1,2), size(temp2,2));
    figure(length(dinfo)+1)
    for seg = 1:numSegs
        flow2_interp = interp1(temp2(:,seg), flow2(:,seg), temp1(:,seg), 'linear', NaN);
        flow_diff = flow1(:,seg) - flow2_interp;
        plot(temp1(:,seg), flow_diff, '-', 'LineWidth', 0.5, 'DisplayName', num2str(seg))
        hold on
    end
    set(gca,'FontSize',15)
    xlabel('Temperature (\circC)','FontSize',20)
    ylabel('\Delta Power (mW)','FontSize',20)
    title(sprintf('Cu2 only', dinfo(2).name, dinfo(1).name), 'Interpreter','none')
    xlim([20 325])
    yline(0, '--k', 'LineWidth', 0.5)  % reference zero line
