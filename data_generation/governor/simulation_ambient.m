clear all
clear global
close % close graphics windows


%% Run TD with pre-set amebint input on mechanical power

td_switch = 0; % if 1 run TD, o.w. 0 just load previous results

dt_const = 0.005;  % uniform time step, o.w. use default in .m sys file.
T_tot = 1200; % total time, if 0 then use default 15s

% n_gen0 = 1; % simulate 1 fault location
% fault_loc = [15 0];
% FO_input = [0.142857 0 0.005];% frequency, phase, amplitude
% % FO_input = [0.5 0 0.005];% frequency, phase, amplitude

dfile = 'data16m_ambient.m'; % pss+damping
% dfile = 'data16m_damp_nopss1.m'; % 2 dampings, no pss, 03/14/23
pathname = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/';

% simu_flag = 0;
if td_switch == 1   
    % Run TD simulation
    simu_flag = s_simu64_ambient(pathname,dfile,dt_const,T_tot);
else
    simu_flag = 1;
end

%% Load simulation results
freq_ambient = cell(1,1);
ang_ambient = cell(1,1);
t_ambient = cell(1,1);
n_gen = zeros(1,1); % number of generators
n_t = zeros(1,1); % number of samples

load_all_data = 1; % if 1, load all raw data, else just load freq and angle

result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/simulation_result/';
if simu_flag == 1
    if load_all_data == 1
        filename = strcat(result_path,'ambient_td_',num2str(T_tot),'s.mat');
%     filename = strcat(result_path,'ambient_td_results.mat');
    else
        filename = strcat(result_path,'ambient_td_',num2str(T_tot),'s_essential.mat');
    end
%     filename = strcat(result_path,'ambient_td_1250s_damp_nopss.mat');
%     filename = strcat(result_path,'ambient_td_1210s_2damp_nopss.mat'); % damp +transient, no pss 
    disp(filename)
    load(filename)
    ang_ambient{1} = mac_ang;
    t_ambient{1} = t;
    n_gen(1) = size(mac_ang,1);
    n_t(1) = size(mac_ang,2);
    freq_ambient{1} = mac_spd - ones(n_gen,n_t); % speed deviation
    disp('Simulated ambient data loaded successfully!')
else
    disp('Load failed!!')
    freq_ambient{1} = 'Nan'; % failed simulation
    ang_ambient{1} = 'Nan';
    t_ambient{1} = 'Nan';
end


%% impz inference using simulated ambient data
% To be finished
% % load ambient_td_results.mat directly
% T_tot = 2000;
% dt = 0.01;
dt = dt_const;
T_start = 20; % we use the data start from 20s

ambient_inference_switch = 0; % if 1, run inference algorithm

if ambient_inference_switch == 1
    disp('Ambient inference for impz start...')
%     [T1,freq_data] = ambient_data_process(freq_ambient{1}(:,1:12000*1),t_ambient{1}(1:12000*1),T_start); % remove first few seconds
    [T1,freq_data] = ambient_data_process(freq_ambient{1},t_ambient{1},T_start); 
    freq_resp = ambient_frequency_response(freq_data',n_gen,dt); % original
%     freq_resp = ambient_frequency_response(freq_data(:,1:12000*15)',n_gen,dt); % 03/09/23, try 1/2min only
else
    disp('Load inferred freq impz instead.')
    filename = strcat(result_path,'ambient_1200s_inferred_freq_resp.mat'); % saved on 01/30/23
    load(filename)
end

% inferred response of input at generator 1
freq_resp1 = freq_resp{1};

disp('Ambient inference for impz finished.')

%% Load impz simulation results
n_gen0 = 16;
simu_flag1 = ones(n_gen0,1);
freq_impz = cell(n_gen0,1);
ang_impz = cell(n_gen0,1);
t_impz = cell(n_gen0,1);
n_gen_impz = zeros(n_gen0,1); % number of generators
n_t_impz = zeros(n_gen0,1); % number of samples

% bus voltage and angle
bus_v_impz = cell(n_gen0,1);
theta_impz = cell(n_gen0,1);

T_tot = 15;
result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/simulation_result/';
for idx = 1 : n_gen0
    mac_fault = idx;
    if simu_flag1(idx) == 1
%         filename = strcat(result_path,'impulse_td_fault',num2str(mac_fault),'.mat');
        filename = strcat(result_path,'impz_td/impulse_td_fault',num2str(mac_fault),'_',num2str(T_tot),'s.mat');
%         filename = strcat(result_path,'impz_td_0309/impulse_td_fault',num2str(mac_fault),'_',num2str(T_tot),'s.mat');
%         filename = strcat(result_path,'impz_td_0310_nopss/impulse_td_fault',num2str(mac_fault),'_',num2str(T_tot),'s.mat');
%         filename = strcat(result_path,'impz_td_0314_damp_nopss/impulse_td_fault',num2str(mac_fault),'_',num2str(T_tot),'s.mat'); % w/ transient damping
        load(filename)
        ang_impz{idx} = mac_ang;
        t_impz{idx} = t;
        n_gen_impz(idx) = size(mac_ang,1);
        n_t_impz(idx) = size(mac_ang,2);
        freq_impz{idx} = mac_spd - ones(n_gen_impz(idx),n_t_impz(idx)); % speed deviation
        % bus voltage and angle
        bus_v_impz{idx} = bus_v;
        theta_impz{idx} = theta;
        if idx == n_gen0
            disp('Simulated impz data loaded successfully.')
        end
    else
        freq_impz{idx} = 'Nan'; % failed simulation
        ang_impz{idx} = 'Nan';
        t_impz{idx} = 'Nan';
    end
end

%% Plot inferred Frequency 

% Time span for prediction
t_range = 5;

T2 = 0 : dt : t_range;
plot_idx = 1 : length(T2);

for idx = 1 : n_gen

freq_resp1 = freq_resp{idx};    
fig1 = figure('DefaultAxesFontSize',18);
freq_resp2 = freq_resp1(plot_idx,:);
move_idx = [9,9];
for i = 1 : n_gen
%     if i ~= 1 && ismember(i,move_idx)
%         freq_resp2(:,i) = -freq_resp2(:,i);
%     end
    freq_resp2(:,i) = freq_resp2(:,i)./max(abs(freq_resp2(:,i)));
    subplot(ceil(sqrt(n_gen)),ceil(sqrt(n_gen)),i)
%     plot(T2,freq_resp2(plot_idx,i),'-','LineWidth',2);
    idx_impz = find(t_impz{1}>1,1)+1;
    temp1 = freq_impz{idx}(i,:)./max(abs(freq_impz{idx}(i,:)));
    plot(t_impz{1}(idx_impz:end)-1,temp1(idx_impz:end),'.-',T2,freq_resp2(plot_idx,i),'-','LineWidth',2);
    xlabel('Time [s]');
    ylabel('scale');
    xlim([0 t_range]);
    title(strcat('\omega ',num2str(i)));
    if i == 1
        legend('simulated','data driven','Location','best');
    end
%     title('Input: \omega_1')
    grid on
end
sgt = sgtitle(strcat('Frequency response: input',num2str(idx)));
sgt.FontSize = 32;
set(fig1,'Position',[10 10 1500 1200])

end

%% PESGM abstract
idx = 1;
freq_resp1 = freq_resp{idx};    
fig1 = figure('DefaultAxesFontSize',18);
freq_resp2 = freq_resp1(plot_idx,:);
gen_idx = [1,7];
for ii = 1 : 2
    nexttile
    i = gen_idx(ii);
    if i ~= 1 && ismember(i,move_idx)
        freq_resp2(:,i) = -freq_resp2(:,i);
    end
    freq_resp2(:,i) = freq_resp2(:,i)./max(abs(freq_resp2(:,i)));
%     subplot(ceil(sqrt(n_gen)),ceil(sqrt(n_gen)),i)
%     plot(T2,freq_resp2(plot_idx,i),'-','LineWidth',2);
    idx_impz = find(t_impz{1}>1,1)+1;
    temp1 = freq_impz{idx}(i,:)./max(abs(freq_impz{idx}(i,:)));
    plot(t_impz{1}(idx_impz:end)-1,temp1(idx_impz:end),'.-',T2,freq_resp2(plot_idx,i),'-','LineWidth',2);
    xlabel('Time [s]');
    ylabel('amplitude');
    xlim([0 t_range]);
    title(strcat('\omega ',num2str(i)));
    if i == 1
        legend('simulated','estimated','Location','best');
    end
%     title('Input: \omega_1')
    grid on
end
% sgt = sgtitle(strcat('Frequency response: input',num2str(idx)));
sgt.FontSize = 32;
set(fig1,'Position',[10 10 600 300])

xx

%% progress review

% idx = 1;
% freq_resp1 = freq_resp{idx};    
% fig1 = figure('DefaultAxesFontSize',18);
% freq_resp2 = freq_resp1(plot_idx,:);
% gen_idx = [1,7];
% for ii = 1 : 9
%     % nexttile
%     i = ii;%gen_idx(ii);
% 
%     idx_impz = find(t_impz{1}>1,1)+1;
%     temp1 = freq_impz{idx}(i,:)./max(abs(freq_impz{idx}(i,:)));
%     plot(t_impz{1}(idx_impz:end)-1,temp1(idx_impz:end),'-','LineWidth',2);
%     hold on
%     xlabel('Time [s]');
%     ylabel('amplitude');
%     xlim([0 15]);
%     grid on
% end
% % sgt = sgtitle(strcat('Frequency response: input',num2str(idx)));
% % sgt.FontSize = 32;
% set(fig1,'Position',[10 10 800 300])


idx = 1;
freq_resp1 = freq_resp{idx};    
fig1 = figure('DefaultAxesFontSize',18);
for ii = 1 : 5
    % nexttile
    i = ii;%gen_idx(ii);
    t00 = 10002;
    plot(t_ambient{1}(t00:end)-500,freq_ambient{1}(i,t00:end),'-','LineWidth',2);
    hold on
    xlabel('Time [s]');
    ylabel('amplitude');
    xlim([0 15]);
    grid on
end
% sgt = sgtitle(strcat('Frequency response: input',num2str(idx)));
% sgt.FontSize = 32;
set(fig1,'Position',[10 10 800 300])




%% FFT test 01/24/2023
figure
fft_diff = cell(n_gen,1);
for idx = 1 : n_gen

test_t1 = t_impz{1}(idx_impz:1200)-1;
test_series1 = freq_impz{1}(idx,:)./max(abs(freq_impz{1}(idx,:)));
test_series1 = test_series1(idx_impz:1200);
test_t2 = T2(3:end);
test_series2 = freq_resp2(plot_idx(3:end),idx);

% figure
% plot(test_t1,test_series1,'.-',test_t2,test_series2,'-');
% title(strcat('impulse response',' Gen',num2str(idx)))
% hold off



Fs = 1/(test_t1(2) - test_t1(1));
L = length(test_t1);      
Y1 = fft(test_series1);
Y2 = fft(test_series2)';

temp_diff = (Y1./L) - (Y2./L);
fft_diff{idx} = temp_diff(1:L/2+1);

% remove small-magnitude transform values 
tol = 1e-2;
Y1(abs(Y1) < tol) = 0;
Y2(abs(Y2) < tol) = 0;

% amplitude
P2_1 = abs(Y1./L);
P1_1 = P2_1(1:L/2+1);
P1_1(2:end-1) = 2.*P1_1(2:end-1);

P2_2 = abs(Y2./L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2.*P1_2(2:end-1);

% angle
PA_11 = angle(Y1);
PA_1 = PA_11(1:L/2+1);
% PA_1(2:end-1) = 2.*PA_1(2:end-1);

PA_22 = angle(Y2);
PA_2 = PA_22(1:L/2+1);
% PA_2(2:end-1) = 2.*PA_2(2:end-1);



f = Fs*(0:(L/2))/L;

% figure
% subplot(2,1,1)
% plot(f,P1_1,'.-',f,P1_2,'-','LineWidth',1.5)
% title(strcat('amplitude',' Gen',num2str(idx)))
% subplot(2,1,2)
% plot(f,PA_1,'.-',f,PA_2,'-','LineWidth',1.5)
% title('angle')
% hold off

subplot(8,4,idx)
plot(f,P1_1,'.-',f,P1_2,'-','LineWidth',1.5)
xlim([0 10]);
title(strcat('amplitude',' Gen',num2str(idx)))

subplot(8,4,16+idx)
plot(f,PA_1,'.-',f,PA_2,'-','LineWidth',1.5)
xlim([0 10]);
title(strcat('angle',' Gen',num2str(idx)))

end

%% plot fft diff
figure
n_spec = size(fft_diff{1},2);
fft_norm = zeros(n_gen,n_spec);
for i = 1 : n_gen
    for j = 1 : n_spec
        fft_norm(i,j) = norm(fft_diff{i}(j));
    end
    plot(fft_norm(i,:),'DisplayName',strcat('Gen',num2str(i)));
%     plot3(f,real(fft_diff{i}),imag(fft_diff{i}),'DisplayName',strcat('Gen',num2str(i)));
    hold on
end
legend('Location','best','NumColumns',2)
xlim([0 10]);
grid on
hold off



%% Plot rotor speed deviation (omega)
figure
for idx = 1 : n_gen0
    subplot(ceil(sqrt(n_gen0)),ceil(sqrt(n_gen0)),idx)
    for i = 1 : n_gen
        plot(t_impz{idx},freq_impz{idx}(i,:),'DisplayName',strcat('Gen',num2str(i)));
        hold on
    end
    xlabel('time [s]')
    ylabel('speed deviation')
    title(strcat('Gen',num2str(idx)))
    if idx == 1
        legend('Location','best','NumColumns',2)
    end
    sgtitle('Rotor speed deviation impz with different input loc')
    hold off
end


% %% FFT analysis
% 
% % set t_start > 20 to remove unstable initial points
% t_start0 = find(t>20,1);
% % t_start0 = 1;
% 
% f_impz = cell(n_gen0,1);
% amp_impz = cell(n_gen0,1);
% 
% % fft_band = [0.1,1];
% for idx = 1 : n_gen0
% %     [f_temp,p_temp] = spectrum_analysis(freq_ambient{idx}(1,:),t_ambient{idx},[]);
%     [f_temp,p_temp] = spectrum_analysis(freq_ambient{idx}(1,t_start0:end),t_ambient{idx}(t_start0:end),[]);
%     f_impz{idx} = zeros(n_gen(idx),length(f_temp));
%     amp_impz{idx} = zeros(n_gen(idx),length(f_temp));
%     for i = 1 : n_gen(idx)
%         [f_impz{idx}(i,:),amp_impz{idx}(i,:)] = spectrum_analysis(freq_ambient{idx}(i,t_start0:end),t_ambient{idx}(t_start0:end),[]);
%     end
% end