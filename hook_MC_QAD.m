function [H_start, H_settle, detail, ind_IL, ind_IR] = hook_MC_QAD(T, Q, FigOn, N, Mode)
% General Function for Measuring Flash DSC Hook. (Quick and Dirty)
% Syntax: [H_start, H_end] = hook_MC(T, Q, FigOn, N, Mode)
% Description: Give it Flash DSC data with a hook. Function will attempt to
% determine whether the hook is at the start or end of the data (heat vs
% cool mode, can also override manually). 
% Inputs: T - Temperature data [C], column vector
%         Q - Heat flow data [mW], column vector
%         FigOn - logical: if true, plots all sorts of figures; default F 
%         N - number of Monte Carlo samples to take, default 1000
%         Mode (optional) - string: "heat" or "cool", if not included, 
%                           will attempt to figure it out automatically
% Outputs: H_start - Hook tip (single value) and noise [val1; val2]
%          H_end - Estimated end of hook [mean; standard deviation] 
%          detail - cell array, any useful info: Mode, LOOCV SP (from dQ),
%                       x_intersect [mean,std]
% Notes on applicability: The hook method assumes that the hook is
% substantial and settles at some semi-linear value. This means that the 
% total heat capacity is high enough, and the heating rate is appropriate 
% to not wash out the hook too much. Also, there should ideally be at least
% 15 degrees between the hook and any large transitions for the code to 
% work. If this is not sufficient, attempts may be made, but modifications 
% to this code could be in order. One more consideration is that some of
% the code assumes it's ok to cut down the LOOCV calulations to a 4th of
% the data. I've unfortuantely found that doesn't work for at least the
% real data noise, so that one starts at 0.67 the data length.

% This version simplifies the pdfs for parameter choices to all be
% rectangular. No LOOCV calculations are used, so the noise values are
% calculated from all the data minus 10 degrees for the hook.

%% Check whether cooling (hook is at right end) or heating (hook is at left end)
% Will need to add nargin syntax when converting to function

% Logic: If I take the average of the Q data, the hook skews the value a
% bit, but the average will still be more representative of the non-hook
% parts of the data. Thus, if both ends of the data are subtracted from the
% calculated average, the larger difference should be on the hook side.
if nargin < 5
    D_avg = mean(Q);
    L_end = abs(Q(1)-D_avg);
    R_end = abs(Q(end)-D_avg);
    if L_end<R_end
        Mode = "cool";
    else
        Mode = "heat";
    end
end

%% Get some random numbers

if nargin < 4
    N = 1000;
end

N_SP = rand(1,N); 
N_SW = rand(1,N); 
N_IP = rand(1,N); 
N_IW = rand(1,N);

%% Check whether figures are on
if nargin < 3
    FigOn = false;
end

%% Switch statement 
switch Mode
%% Case: Heat   
   case 'heat'
%% PDF and CDF for Settling Point (SP) + Sampling 
% Rules: Calculate the derivative, its noise fit. Find the first value from 
% the hook direction that crosses the noise fit. That's the far point I'll 
% take, and then I'll sample another 5 degrees from that. 
% The noise fit is the fit to all the data after 10 degrees from the hook.

% Get a derivative to work with:
dT = zeros(length(T)-1,1);
dQ = zeros(length(T)-1,1);
for j = 1:length(T)-1
    dT(j) = (T(j+1)+T(j))/2;
    dQ(j) = (Q(j+1)-Q(j))/(T(j+1)-T(j));
end

ind_noise = dT>T(1)+10;
dT_noise = dT(ind_noise);
P_noise = polyfit(dT_noise,dQ(ind_noise),2);

dQ_noise = polyval(P_noise,dT_noise);
RMSE_dNoise = sqrt(sum((dQ(ind_noise)-dQ_noise).^2)/length(dT_noise));

% Inflection point is nice to have for intermediate steps in getting PDF
[~,infl] = min(dQ);   % an index
T_infl = dT(infl);
% Find the first point where dQ crosses over the noise fit:
settles = find(dQ-polyval(P_noise,dT) > 0 & dT > T_infl); % want first val
T_SP_min = T(settles(1));

% Just take the next 5 degrees from there:
T_SP_max = T_SP_min+5;

% Define rectangular pdf:
pdf_SP = zeros(size(T));
ind_SP = T < T_SP_max & T > T_SP_min;
pdf_SP(ind_SP) = 1/(T_SP_max-T_SP_min);

% Sample values:
v_SP = N_SP*(T_SP_max-T_SP_min)+T_SP_min;

% Get closest points to samples values
ind_SP = zeros(size(N_SP));
T_SP = zeros(size(N_SP));
for i = 1:(length(N_SP))
    [~,ind_SP(i)] = min(abs(T-v_SP(i)));
    T_SP(i) = T(ind_SP(i));
end

%% PDF and CDF for Settling Region Width (SW)
% Defined as between 5-10 degrees

w_min = 5;
w_max = 10;

% Really easy analytical solution and inversion
pdf_SW = 1/(w_max-w_min);
% cdf_SW = (w_SW-w_min)/(w_max-w_min);

% Get the sampled widths
w_SW = N_SW*(w_max-w_min)+w_min;

% Going to take it one step farther to get the index of the end point
% defined by T_SP + w_SW . Here's where the closest index matters.

ind_SW = zeros(size(N_SW));
T_SW = zeros(size(N_SW));
for i = 1:(length(N_SW))
    [~,ind_SW(i)] = min(abs(T-(T_SP(i) + w_SW(i))));
    T_SW(i) = T(ind_SW(i));
end

%% PDF and CDF for Inflection Point (IP)
% Rules: From the inflection point, take the 2sigma around it based on the 
% previous definition of sigma.  

infls = dQ<(min(dQ)+2*RMSE_dNoise);  % Only actually ~7 points...
s_IP = max(abs(dT(infls)-T_infl))/4;   % 1 sigma is too much, this can be 1/4th of a sigma.

T_IP_min = T_infl-2*s_IP;
T_IP_max = T_infl+2*s_IP;

% Define rectangular pdf:
pdf_IP = zeros(size(T));
ind_IP = T < T_IP_max & T > T_IP_min;
pdf_IP(ind_IP) = 1/(T_IP_max-T_IP_min);

% Sample values:

v_IP = N_IP.*(T_IP_max-T_IP_min)+T_IP_min;
disp('v_IP'); disp(size(v_IP))
% Get closest points to samples values

ind_IP = zeros(size(N_IP));
T_IP = zeros(size(N_IP));
for i = 1:(length(N_IP))
    disp('v_IP'); disp(v_IP)

    [~,ind_IP(i)] = min(abs(T-v_IP(i)));
    T_IP(i) = T(ind_IP(i));
end

%% PDF and CDF for Inflection Region (Half) Width (IW)
% Rules: Min is the half the range of IP (2s_IP) 
% Max is the smaller distance from the inflection point to the 1/4 height 
% crossover point  

% Get the base of the inflection peak
y_base = polyval(P_noise,dT);
% Get the 1/4 peak range
HM = dQ(infl)+(y_base(infl)-dQ(infl))*3/4;
get_HM = find(dQ<HM);  
ws = [abs(min(dT(get_HM))-T_infl),abs(max(dT(get_HM))-T_infl)];
w_IW_max = min(ws);

% Get the minimum width
w_IW_min = 2*s_IP;

% Really easy analytical solution and inversion
pdf_IW = 1/(w_IW_max-w_IW_min);

% Get the sampled widths
w_IW = N_IW*(w_IW_max-w_IW_min)+w_IW_min;

% As a last step, let's get the indexes around the sampled inflection pts
ind_IL = zeros(size(T));  % Lower T than the inflection point
ind_IR = zeros(size(T));  % Higher T than the inflection point
T_ILl = zeros(size(T));
T_IRr = zeros(size(T));

for i = 1:length(N_IW)
    [~,ind_IL(i)] = min(abs(T-(T_IP(i)-w_IW(i)))); % Closest T index to Tinfl-width
    [~,ind_IR(i)] = min(abs(T-(T_IP(i)+w_IW(i)))); % Closest T index to Tinfl+width
    T_ILl(i) = T(ind_IL(i));
    T_IRr(i) = T(ind_IR(i));
end

%% So now to get the intersections.
% Will do the calculations here for now, since I have the indexes right
% here, but possibly think about making a separate function later on.

y_intersects = zeros(N,1);
x_intersects = zeros(N,1);
for i = 1:N
    % Settled region (SR)
    T_SR = T(ind_SP(i):ind_SW(i));
    Q_SR = Q(ind_SP(i):ind_SW(i));
    P_SR = polyfit(T_SR,Q_SR,1);
    m_SR = P_SR(1);
    b_SR = P_SR(2);
    
    % Inflection region (IR)
    if abs(ind_IL(i) - ind_IR(i)) < 1
        x2 = T(ind_IL(i)+1);
        x1 = T(ind_IL(i));
        y2 = Q(ind_IL(i)+1);
        y1 = Q(ind_IL(i));
        T_IR = [x1 x2];
        Q_IR = [y1 y2];
        P_IR = polyfit(T_IR,Q_IR,1);
        m_IR = P_IR(1);
        b_IR = P_IR(2);

        x_intersects(i) = (b_IR-b_SR)/(m_SR-m_IR);
        y_intersects(i) = m_SR*x_intersects(i)+b_SR;
    else
    T_IR = T(ind_IL(i):ind_IR(i));
    Q_IR = Q(ind_IL(i):ind_IR(i));
    P_IR = polyfit(T_IR,Q_IR,1);
    m_IR = P_IR(1);
    b_IR = P_IR(2);
    
    x_intersects(i) = (b_IR-b_SR)/(m_SR-m_IR);
    y_intersects(i) = m_SR*x_intersects(i)+b_SR;
    end
end

%% Case: Cool  
    case 'cool'
%% PDF and CDF for Settling Point (SP) + Sampling 

% Get a derivative to work with:
dT = zeros(length(T)-1,1);
dQ = zeros(length(T)-1,1);
for j = 1:length(T)-1
    dT(j) = (T(j+1)+T(j))/2;
    dQ(j) = (Q(j+1)-Q(j))/(T(j+1)-T(j));
end

ind_noise = dT<T(end)-10;
dT_noise = dT(ind_noise);
P_noise = polyfit(dT_noise,dQ(ind_noise),2);
dQ_noise = polyval(P_noise,dT_noise);
RMSE_dNoise = sqrt(sum((dQ(ind_noise)-dQ_noise).^2)/length(dT_noise));

% Inflection point is nice to have for intermediate steps in getting PDF
[~,infl] = min(dQ);   % an index
T_infl = dT(infl);

% Find the first point where dQ crosses over the noise fit:
settles = find(dQ-polyval(P_noise,dT) > 0 & dT < T_infl); % want last val
T_SP_max = T(settles(end));

% Just take the next 5 degrees from there:
T_SP_min = T_SP_max-5;

% Define rectangular pdf:
pdf_SP = zeros(size(T));
ind_SP = T < T_SP_max & T > T_SP_min;
pdf_SP(ind_SP) = 1/(T_SP_max-T_SP_min);

% Sample values:
v_SP = N_SP*(T_SP_max-T_SP_min)+T_SP_min;

% Get closest points to samples values
ind_SP = zeros(size(N_SP));
T_SP = zeros(size(N_SP));
for i = 1:(length(N_SP))
    [~,ind_SP(i)] = min(abs(T-v_SP(i)));
    T_SP(i) = T(ind_SP(i));
end

%% PDF and CDF for Settling Region Width (SW)

% Actual temperature chosen is dependent on the sampled T_SP
% But I'm going to sample the uniform distributions of widths

w_min = 5;
w_max = 10;

% Really easy analytical solution and inversion
pdf_SW = 1/(w_max-w_min);
% cdf_SW = (w_SW-w_min)/(w_max-w_min);

% Get the sampled widths
w_SW = N_SW*(w_max-w_min)+w_min;

% Going to take it one step farther to get the index of the end point
% defined by T_SP + w_SW . Here's where the closest index matters.

ind_SW = zeros(size(N_SW));
T_SW = zeros(size(N_SW));
for i = 1:(length(N_SW))
    [~,ind_SW(i)] = min(abs(T-(T_SP(i) - w_SW(i))));
    T_SW(i) = T(ind_SW(i));
end

%% PDF and CDF for Inflection Point (IP)

infls = dQ<(min(dQ)+2*RMSE_dNoise);  % Only actually ~7 points...
s_IP = max(abs(dT(infls)-T_infl))/4;   % 1 sigma is too much, this can be 1/4th of a sigma.

T_IP_min = T_infl-2*s_IP;
T_IP_max = T_infl+2*s_IP;

% Define rectangular pdf:
pdf_IP = zeros(size(T));
ind_IP = T < T_IP_max & T > T_IP_min;
pdf_IP(ind_IP) = 1/(T_IP_max-T_IP_min);

% Sample values:
v_IP = N_IP*(T_IP_max-T_IP_min)+T_IP_min;

% Get closest points to samples values
ind_IP = zeros(size(N_IP));
T_IP = zeros(size(N_IP));
for i = 1:(length(N_IP))
    [~,ind_IP(i)] = min(abs(T-v_IP(i)));
    T_IP(i) = T(ind_IP(i));
end

%% PDF and CDF for Inflection Region (Half) Width (IW)

% Get the base of the inflection peak
y_base = polyval(P_noise,dT);
% Get the 1/4 peak range
HM = dQ(infl)+(y_base(infl)-dQ(infl))*3/4;
get_HM = find(dQ<HM);  
ws = [abs(min(dT(get_HM))-T_infl),abs(max(dT(get_HM))-T_infl)];
w_IW_max = min(ws);

% Get the minimum width
w_IW_min = 2*s_IP;

% Really easy analytical solution and inversion
pdf_IW = 1/(w_IW_max-w_IW_min);

% Get the sampled widths
w_IW = N_IW*(w_IW_max-w_IW_min)+w_IW_min;

% As a last step, let's get the indexes around the sampled inflection pts
ind_IL = zeros(size(T));  % Lower T than the inflection point
ind_IR = zeros(size(T));  % Higher T than the inflection point
T_ILl = zeros(size(T));
T_IRr = zeros(size(T));

for i = 1:length(N_IW)
    [~,ind_IL(i)] = min(abs(T-(T_IP(i)-w_IW(i)))); % Closest T index to Tinfl-width
    [~,ind_IR(i)] = min(abs(T-(T_IP(i)+w_IW(i)))); % Closest T index to Tinfl+width
    T_ILl(i) = T(ind_IL(i));
    T_IRr(i) = T(ind_IR(i));
end

%% So now to get the intersections.
% Will do the calculations here for now, since I have the indexes right
% here, but possibly think about making a separate function later on.

y_intersects = zeros(N,1);
x_intersects = zeros(N,1);
for i = 1:N
    % Settled region (SR)
    T_SR = T(ind_SW(i):ind_SP(i));
    Q_SR = Q(ind_SW(i):ind_SP(i));
    P_SR = polyfit(T_SR,Q_SR,1);
    m_SR = P_SR(1);
    b_SR = P_SR(2);
    
    % Inflection region (IR)
    T_IR = T(ind_IL(i):ind_IR(i));
    Q_IR = Q(ind_IL(i):ind_IR(i));
    P_IR = polyfit(T_IR,Q_IR,1);
    m_IR = P_IR(1);
    b_IR = P_IR(2);
    
    x_intersects(i) = (b_IR-b_SR)/(m_SR-m_IR);
    y_intersects(i) = m_SR*x_intersects(i)+b_SR;
end
       
%% End Switch Statement
end

%% Define Outputs: mean and standard deviation of intersects

% Actual mean and std dev of the y-intersects:
m_Yint = mean(y_intersects);
s_Yint = std(y_intersects);

%% Define Outputs: end of data and associated RMSE noise estimate
switch Mode
    case 'heat'
    % The noise fit is the fit to all the data after 10 degrees from the hook.
    ind_noise = T>T(1)+10;
    T_noise = T(ind_noise);
    Pr = polyfit(T_noise,Q(ind_noise),3);
    Q_noise = polyval(Pr,T_noise);
    RMSE_noise = sqrt(sum((Q(ind_noise)-Q_noise).^2)/length(T_noise));
    bottom = Q(1);
    % if FigOn
    %     figure;
    %     hold on
    %     plot(T,Q)
    %     plot(T_noise,Q_noise)
    %     xlabel('Temperature (C)')
    %     ylabel('Heat Flow (mW)')
    %     title('Cubic Noise Fit')
    % end  
    
    case 'cool'
    ind_noise = T<T(end)-10;
    T_noise = T(ind_noise);
    Pr = polyfit(T_noise,Q(ind_noise),3);
    Q_noise = polyval(Pr,T_noise);
    RMSE_noise = sqrt(sum((Q(ind_noise)-Q_noise).^2)/length(T_noise));
    bottom = Q(end);
    % if FigOn
    %     figure;
    %     hold on
    %     plot(T,Q)
    %     plot(T_noise,Q_noise)
    %     xlabel('Temperature (C)')
    %     ylabel('Heat Flow (mW)')
    %     title('Cubic Noise Fit')
    % end  
end

%% Save to output variables

H_start = [bottom; RMSE_noise];
H_settle = [m_Yint; s_Yint];
detail = {Mode, dT_noise,[mean(x_intersects),std(x_intersects)]};

%% Plotting figures:

if FigOn
    switch Mode
        case 'heat'
        %% Settling Point Composite Figure
%         limX = [min(T),max(T)];
%         limY1 = [min(Q), max(Q)];
%         limY2 = [min(dQ), max(dQ)];
%         limY3 = [0,max(pdf_SP)*1.1];
% 
%         figure;
%         subplot(3,1,1)
%         hold on
%         plot(T,Q)
%         plot(T_SP_min*ones(1,2),limY1,'k--')
%         plot(T_SP_max*ones(1,2),limY1,'k--')
% 
%         ylabel('Heat Flow (mW)')
%         xlim(limX)
%         ylim(limY1)
% 
%         subplot(3,1,2)
%         hold on
%         plot(dT,dQ)
%         plot(T_SP_min*ones(1,2),limY2,'k--')
%         plot(T_SP_max*ones(1,2),limY2,'k--')
%         plot(dT,polyval(P_noise,dT))
%         plot(dT,polyval(P_noise,dT)+2*RMSE_dNoise,'c--')
%         plot(dT,polyval(P_noise,dT)-2*RMSE_dNoise,'c--')
% 
%         ylabel('dQ/dT')
%         xlim(limX)
%         ylim(limY2)
%         h = zeros(1,1);
%         h(1) = plot(NaN,NaN,'c--');
%         legend(h,'+/- 2 RMSE from noise fit')
% 
%         subplot(3,1,3)
%         hold on
%         plot(T,pdf_SP)
%         plot(T_SP_min*ones(1,2),limY3,'k--')
%         plot(T_SP_max*ones(1,2),limY3,'k--')
% 
%         xlim(limX)
%         ylim(limY3)
%         ylabel('Settling point PDF')
%         xlabel('Temperature (C)')
% 
%         %% Settling width composite figure     
%         Tmean = mean(T_SP);
% 
%         figure;
%         subplot(3,1,1)
%         hold on
%         plot(T,Q)
%         plot(Tmean*ones(1,2),limY1,'g--')
%         plot((Tmean+w_min)*ones(1,2),limY1,'k--')
%         plot((Tmean+w_max)*ones(1,2),limY1,'k--')
% 
%         ylabel('Heat Flow (mW)')
%         xlim(limX)
%         ylim(limY1)
% 
%         subplot(3,1,2)
%         hold on
%         plot(dT,dQ)
%         plot(dT,polyval(P_noise,dT))
%         plot(Tmean*ones(1,2),limY2,'g--')
%         plot((Tmean+w_min)*ones(1,2),limY2,'k--')
%         plot((Tmean+w_max)*ones(1,2),limY2,'k--')
% 
%         ylabel('dQ/dT')
%         xlim(limX)
%         ylim(limY2)
% 
%         subplot(3,1,3)
%         hold on
%         plot(T,pdf_SP,'g--')
%         plot(Tmean+[w_min,w_max],pdf_SW*ones(1,2),'b')
%         plot(Tmean*ones(1,2),limY3,'g--')
%         plot((Tmean+w_min)*ones(1,2),limY3,'k--')
%         plot((Tmean+w_max)*ones(1,2),limY3,'k--')
% 
%         xlim(limX)
%         ylim(limY3)
%         ylabel('Settling width PDF')
%         xlabel('Temperature (C)')   
% 
%         %% Inflection point composite figure
% 
%         limX = [min(T),max(T)];
%         limY1 = [min(Q), max(Q)];
%         limY2 = [min(dQ), max(dQ)];
% %         limY3 = [0,max(pdf_IP)*1.1]; 
%         limY3 = [0 1.1];
%         figure;
%         subplot(3,1,1)
%         hold on
%         plot(T,Q)
%         plot(T_IP_min*ones(1,2),limY1,'k--')
%         plot(T_IP_max*ones(1,2),limY1,'k--')
% 
%         ylabel('Heat Flow (mW)')
%         xlim(limX)
%         ylim(limY1)
% 
%         subplot(3,1,2)
%         hold on
%         plot(dT,dQ)
%         plot(T_IP_min*ones(1,2),limY2,'k--')
%         plot(T_IP_max*ones(1,2),limY2,'k--')
% 
%         ylabel('dQ/dT')
%         xlim(limX)
%         ylim(limY2)
% 
%         subplot(3,1,3)
%         hold on
%         plot(T,pdf_IP)
%         plot(T_IP_min*ones(1,2),limY3,'k--')
%         plot(T_IP_max*ones(1,2),limY3,'k--')
% 
%         xlim(limX)
%         ylim(limY3)
%         ylabel('Inflection point PDF')
%         xlabel('Temperature (C)')
% 
%         %% Inflection region composite figure
%         m_IP = mean(T_IP);
% 
%         figure;
%         subplot(3,1,1)
%         hold on
%         plot(T,Q)
%         plot(T_IP_min*ones(1,2),limY1,'g--')
%         plot(T_IP_max*ones(1,2),limY1,'g--')
%         plot((m_IP-w_IW_min)*ones(1,2),limY1,'k--')
%         plot((m_IP-w_IW_max)*ones(1,2),limY1,'k--')
%         plot((m_IP+w_IW_min)*ones(1,2),limY1,'k:')
%         plot((m_IP+w_IW_max)*ones(1,2),limY1,'k:')
% 
%         ylabel('Heat Flow (mW)')
%         xlim(limX)
%         ylim(limY1)
% 
%         subplot(3,1,2)
%         hold on
%         plot(dT,dQ)
%         plot(dT,HM*ones(size(dT)),'b:')
%         plot(T_IP_min*ones(1,2),limY2,'g--')
%         plot(T_IP_max*ones(1,2),limY2,'g--')
%         plot((m_IP-w_IW_min)*ones(1,2),limY2,'k--')
%         plot((m_IP-w_IW_max)*ones(1,2),limY2,'k--')
%         plot((m_IP+w_IW_min)*ones(1,2),limY2,'k:')
%         plot((m_IP+w_IW_max)*ones(1,2),limY2,'k:')        
% 
%         ylabel('dQ/dT')
%         xlim(limX)
%         ylim(limY2)
% 
%         subplot(3,1,3)
%         hold on
%         plot(T,pdf_IP,'g')
%         plot(T_IP_min*ones(1,2),limY3,'g--')
%         plot(T_IP_max*ones(1,2),limY3,'g--')
%         plot((m_IP-w_IW_min)*ones(1,2),limY3,'k--')
%         plot((m_IP-w_IW_max)*ones(1,2),limY3,'k--')
%         plot((m_IP+w_IW_min)*ones(1,2),limY3,'k:')
%         plot((m_IP+w_IW_max)*ones(1,2),limY3,'k:')       
%         plot(m_IP-[w_IW_min,w_IW_max],pdf_IW*ones(1,2),'b')
%         plot(m_IP+[w_IW_min,w_IW_max],pdf_IW*ones(1,2),'b--')
% 
%         xlim(limX)
%         ylim(limY3)
%         ylabel('Inflection region PDF')
%         xlabel('Temperature (C)')

        %% Chosen values with comparison to previous method
        
        hook_T = T(1);  
        hook_var(T, Q, hook_T, 0.75, 10, 5, true, true);
        % ^ Makes a figure with construction lines
        % hold on
        % plot(x_intersects,y_intersects,'o','MarkerEdgeColor',[128/255 128/255 128/255])
        % plot(T,Q,'k')
        % plot(T(ind_SW),Q(ind_SW),'go')
        % plot(T(ind_SP),Q(ind_SP),'mo')
        % plot(T(ind_IL),Q(ind_IL),'bo')
        % plot(T(ind_IR),Q(ind_IR),'bo')
        % plot(T(ind_IP),Q(ind_IP),'ro')
        % xlabel('Temperature (C)')
        % ylabel('Heat Flow (mW)')
     
        case 'cool'
        %% Settling Point composite figure
        
        limX = [min(T),max(T)];
        limY1 = [min(Q), max(Q)];
        limY2 = [min(dQ), max(dQ)];
        limY3 = [0,max(pdf_SP)*1.1];
        
        % figure;
        % subplot(3,1,1)
        % hold on
        % plot(T,Q)
        % plot(T_SP_min*ones(1,2),limY1,'k--')
        % plot(T_SP_max*ones(1,2),limY1,'k--')
        % 
        % ylabel('Heat Flow (mW)')
        % xlim(limX)
        % ylim(limY1)
        % 
        % subplot(3,1,2)
        % hold on
        % plot(dT,dQ)
        % plot(T_SP_min*ones(1,2),limY2,'k--')
        % plot(T_SP_max*ones(1,2),limY2,'k--')
        % plot(dT,polyval(P_noise,dT))
        % plot(dT,polyval(P_noise,dT)+2*RMSE_dNoise,'c--')
        % plot(dT,polyval(P_noise,dT)-2*RMSE_dNoise,'c--')
        % 
        % ylabel('dQ/dT')
        % xlim(limX)
        % ylim(limY2)
        % h = zeros(1,1);
        % h(1) = plot(NaN,NaN,'c--');
        % legend(h,'+/- 2 RMSE from noise fit')
        % 
        % subplot(3,1,3)
        % hold on
        % plot(T,pdf_SP)
        % plot(T_SP_min*ones(1,2),limY3,'k--')
        % plot(T_SP_max*ones(1,2),limY3,'k--')
        % 
        % xlim(limX)
        % ylim(limY3)
        % ylabel('Settling point PDF')
        % xlabel('Temperature (C)')
        % 
        % %% Settling width composite figure
        % 
        % Tmean = mean(T_SP);
        % 
        % figure;
        % subplot(3,1,1)
        % hold on
        % plot(T,Q)
        % plot(Tmean*ones(1,2),limY1,'g--')
        % plot((Tmean-w_min)*ones(1,2),limY1,'k--')
        % plot((Tmean-w_max)*ones(1,2),limY1,'k--')
        % 
        % ylabel('Heat Flow (mW)')
        % xlim(limX)
        % ylim(limY1)
        % 
        % subplot(3,1,2)
        % hold on
        % plot(dT,dQ)
        % plot(dT,polyval(P_noise,dT))
        % plot(Tmean*ones(1,2),limY2,'g--')
        % plot((Tmean-w_min)*ones(1,2),limY2,'k--')
        % plot((Tmean-w_max)*ones(1,2),limY2,'k--')
        % 
        % ylabel('dQ/dT')
        % xlim(limX)
        % ylim(limY2)
        % 
        % subplot(3,1,3)
        % hold on
        % plot(T,pdf_SP,'g--')
        % plot(Tmean+[w_min,w_max],pdf_SW*ones(1,2),'b')
        % plot(Tmean*ones(1,2),limY3,'g--')
        % plot((Tmean-w_min)*ones(1,2),limY3,'k--')
        % plot((Tmean-w_max)*ones(1,2),limY3,'k--')
        % 
        % xlim(limX)
        % ylim(limY3)
        % ylabel('Settling width PDF')
        % xlabel('Temperature (C)')   
        % 
        % %% Inflection point composite figure
        % 
        % limX = [min(T),max(T)];
        % limY1 = [min(Q), max(Q)];
        % limY2 = [min(dQ), max(dQ)];
        % limY3 = [0,max(pdf_IP)*1.1];
        % 
        % figure;
        % subplot(3,1,1)
        % hold on
        % plot(T,Q)
        % plot(T_IP_min*ones(1,2),limY1,'k--')
        % plot(T_IP_max*ones(1,2),limY1,'k--')
        % 
        % ylabel('Heat Flow (mW)')
        % xlim(limX)
        % ylim(limY1)
        % 
        % subplot(3,1,2)
        % hold on
        % plot(dT,dQ)
        % plot(T_IP_min*ones(1,2),limY2,'k--')
        % plot(T_IP_max*ones(1,2),limY2,'k--')
        % 
        % ylabel('dQ/dT')
        % xlim(limX)
        % ylim(limY2)
        % 
        % subplot(3,1,3)
        % hold on
        % plot(T,pdf_IP)
        % plot(T_IP_min*ones(1,2),limY3,'k--')
        % plot(T_IP_max*ones(1,2),limY3,'k--')
        % 
        % xlim(limX)
        % ylim(limY3)
        % ylabel('Inflection point PDF')
        % xlabel('Temperature (C)')
        % 
        % %% Inflection region composite figure
        % 
        %  m_IP = mean(T_IP);
        % 
        % figure;
        % subplot(3,1,1)
        % hold on
        % plot(T,Q)
        % plot(T_IP_min*ones(1,2),limY1,'g--')
        % plot(T_IP_max*ones(1,2),limY1,'g--')
        % plot((m_IP+w_IW_min)*ones(1,2),limY1,'k--')
        % plot((m_IP+w_IW_max)*ones(1,2),limY1,'k--')
        % plot((m_IP-w_IW_min)*ones(1,2),limY1,'k:')
        % plot((m_IP-w_IW_max)*ones(1,2),limY1,'k:')
        % 
        % ylabel('Heat Flow (mW)')
        % xlim(limX)
        % ylim(limY1)
        % 
        % subplot(3,1,2)
        % hold on
        % plot(dT,dQ)
        % plot(dT,HM*ones(size(dT)),'b:')
        % plot(T_IP_min*ones(1,2),limY2,'g--')
        % plot(T_IP_max*ones(1,2),limY2,'g--')
        % plot((m_IP+w_IW_min)*ones(1,2),limY2,'k--')
        % plot((m_IP+w_IW_max)*ones(1,2),limY2,'k--')
        % plot((m_IP-w_IW_min)*ones(1,2),limY2,'k:')
        % plot((m_IP-w_IW_max)*ones(1,2),limY2,'k:')        
        % 
        % ylabel('dQ/dT')
        % xlim(limX)
        % ylim(limY2)
        % 
        % subplot(3,1,3)
        % hold on
        % plot(T,pdf_IP,'g')
        % plot(T_IP_min*ones(1,2),limY3,'g--')
        % plot(T_IP_max*ones(1,2),limY3,'g--')
        % plot((m_IP+w_IW_min)*ones(1,2),limY3,'k--')
        % plot((m_IP+w_IW_max)*ones(1,2),limY3,'k--')
        % plot((m_IP-w_IW_min)*ones(1,2),limY3,'k:')
        % plot((m_IP-w_IW_max)*ones(1,2),limY3,'k:')       
        % plot(m_IP+[w_IW_min,w_IW_max],pdf_IW*ones(1,2),'b')
        % plot(m_IP-[w_IW_min,w_IW_max],pdf_IW*ones(1,2),'b--')
        % 
        % xlim(limX)
        % ylim(limY3)
        % ylabel('Inflection region PDF')
        % xlabel('Temperature (C)')
        
        %% Chosen values with comparison to previous method
        
        hook_T = T(end);        
        hook_var(T, Q, hook_T, 0.75, 10, 5, false, true);
        % ^ Makes a figure with construction lines
        hold on
        % plot(x_intersects,y_intersects,'o','MarkerEdgeColor',[128/255 128/255 128/255])
        % plot(T,Q,'k')
        % plot(T(ind_SW),Q(ind_SW),'go')
        % plot(T(ind_SP),Q(ind_SP),'mo')
        % plot(T(ind_IL),Q(ind_IL),'bo')
        % plot(T(ind_IR),Q(ind_IR),'bo')
        % plot(T(ind_IP),Q(ind_IP),'ro')
        % xlabel('Temperature (C)')
        % ylabel('Heat Flow (mW)')
        
    end
    
    %% Histogran of y intersects and variable effects
    % figure;
    % hold on
    % h = histogram(y_intersects);
    % title('y-intersects')
    % y_h = (h.Values)';
    % x_h = (h.BinEdges(1:end-1)+h.BinWidth/2)';
    % F_h = fit(x_h, y_h, 'gauss1');
    % plot(x_h,F_h(x_h),'r')
    % 
    % figure;
    % subplot(2,2,1)
    % plot(T_SP,y_intersects,'o')
    % xlabel('Settling Point')
    % ylabel('Y intersect')
    % 
    % subplot(2,2,2)
    % plot(w_SW,y_intersects,'o')
    % xlabel('Settling Width')
    % ylabel('Y intersect')
    % 
    % subplot(2,2,3)
    % plot(T_IP,y_intersects,'o')
    % xlabel('Inflection Point')
    % ylabel('Y intersect')
    % 
    % subplot(2,2,4)
    % plot(w_IW,y_intersects,'o')
    % xlabel('Inflection Width')
    % ylabel('Y intersect')

end

end
