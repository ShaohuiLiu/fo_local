clear all
clear global
close % close graphics windows

% n_gen0 = 16;
n_exc = 9;
n_gen0 = n_exc;

% simu_type = 'impz'; % V_ref
% simu_type = 'xxxx'; % V_R
simu_type = 'xxxy'; % E_fd

%% Run TD for impulse on all exciters

T_tot = 15;

td_switch = 1; % if 1 run TD, o.w. 0 just load previous results

dt_const = 0.005;  % uniform time step, o.w. use default in .m sys file. 

% dfile = 'data16m.m';
% saved separate file for exciter simulation
dfile = 'data16m_ambient_v.m'; % additional damping added simulated on 03/08/23 
% dfile = 'data16m_damp_nopss.m'; % additional damping no PSS simulated on 03/10/23
% dfile = 'data16m_damp_nopss1.m'; % additional damping w/ transient damping no PSS simulated on 03/14/23
pathname = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests_V/';

simu_flag = zeros(n_gen0,1);
if td_switch == 1
    for i = 1 : n_exc
        mac_fault = i;% fault location    
        % Run TD simulation
        simu_flag(i) = s_simu64_exciter(pathname,dfile,mac_fault,dt_const,T_tot,simu_type); % run 15s
    end
else
    simu_flag = ones(n_gen0,1);
end


%% Load simulation results
freq_impz = cell(n_gen0,1);
ang_impz = cell(n_gen0,1);
theta_impz = cell(n_gen0,1);
bus_v_impz = cell(n_gen0,1);
t_impz = cell(n_gen0,1);
n_gen = zeros(n_gen0,1); % number of generators
n_t = zeros(n_gen0,1); % number of samples

% gen voltage vals
ed_impz = cell(n_gen0,1);
eq_impz = cell(n_gen0,1);
edprime_impz = cell(n_gen0,1);
eqprime_impz = cell(n_gen0,1);
Efd_impz = cell(n_gen0,1);


% load earlier results
% result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/simulation_result/impz_td/';
% result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/simulation_result/impz_td_0309/';
% result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/simulation_result/impz_td_0310_nopss/';
result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests_V/simulation_result/impz_td_vr_0816/';
result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests_V/simulation_result/impz_td_Efd2_0823/';
for idx = 1 : n_gen0
    mac_fault = idx;
    if simu_flag(idx) == 1
        % filename = strcat('impulse_td_exciter_fault',num2str(mac_fault),'_',num2str(T_tot),'s.mat');
        % filename = strcat(result_path,'impulse_td_fault',num2str(mac_fault),'.mat');
        % filename = strcat(result_path,'impulse_td_exciter_VR_fault',num2str(mac_fault),'_',num2str(T_tot),'s.mat');
        filename = strcat(result_path,'impulse_td_exciter_Efd_fault',num2str(mac_fault),'_',num2str(T_tot),'s.mat');
        load(filename)
        ang_impz{idx} = mac_ang;
        t_impz{idx} = t;
        n_gen(idx) = size(mac_ang,1);
        n_t(idx) = size(mac_ang,2);
        freq_impz{idx} = mac_spd - ones(n_gen(idx),n_t(idx)); % speed deviation
        bus_v_impz{idx} = bus_v;
        theta_impz{idx} = theta;
        % other voltage vals
        ed_impz{idx} = ed;
        eq_impz{idx} = eq;
        edprime_impz{idx} = edprime;
        eqprime_impz{idx} = eqprime;
        Efd_impz{idx} = Efd;
    else
        disp('No read!')
        freq_impz{idx} = 'Nan'; % failed simulation
        ang_impz{idx} = 'Nan';
        t_impz{idx} = 'Nan';
    end
    fprintf('Sys %d loaded. \n',idx)
end

xx

%% FFT analysis
f_impz = cell(n_gen0,1);
amp_impz = cell(n_gen0,1);

fft_band = [0.1,0.4];
for idx = 1 : n_gen0
    disp('Determine size.')
    [f_temp,p_temp] = spectrum_analysis(freq_impz{idx}(1,:),t_impz{idx},[]); % check size
    f_impz{idx} = zeros(n_gen(idx),length(f_temp));
    amp_impz{idx} = zeros(n_gen(idx),length(f_temp));
    fprintf('Analyze input %d with band [%.1f,%.1f]. \n',idx,fft_band(1),fft_band(2));
    for i = 1 : n_gen(idx)
        % Use default band [0.1,Fs]
        [f_impz{idx}(i,:),amp_impz{idx}(i,:)] = spectrum_analysis(freq_impz{idx}(i,:),t_impz{idx},[]);
%         % Use designated band
%         [f_impz{idx}(i,:),amp_impz{idx}(i,:)] = spectrum_analysis(freq_impz{idx}(i,:),t_impz{idx},fft_band);
    end
end

%% plot spectrum
plt_switch = 1; % if 1 plot

if plt_switch == 1
    % fft of rotor speed deviation (omega)
    figure
    for idx = 1 : n_gen0
        subplot(4,4,idx)
        for i = 1 : n_gen
            if i == idx
                plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'LineWidth',2,'DisplayName',strcat('Gen',num2str(i)));
                hold on
            else
                plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'DisplayName',strcat('Gen',num2str(i)));
                hold on
            end
        end
        xlabel('frequency [Hz]')
        ylabel('amplitude')
%         xlim([-0.5,10])
        xlim([-0.1,3])
        title(strcat('Gen',num2str(idx)))
        if idx == 1
            legend('Location','best','NumColumns',2)
        end
        sgtitle('FFT of \omega impz with different input loc''s')
        hold off
    end

else
    disp('No plot of spectrum.')
end


%% plot spectrum summation to identify modes

amp_impz_sum = zeros(size(amp_impz{1}(1,:)));
for idx = 1 : n_gen0
    for i = 1 : n_gen
        amp_impz_sum = amp_impz_sum + amp_impz{idx}(i,:);
    end
end
figure
plot(f_impz{1}(1,:),amp_impz_sum,'LineWidth',2)
xlabel('frequency [Hz]')
ylabel('summed amplitude')
xlim([-0.1,3])
title('Aggregated spectrum')
hold off

figure
for idx = 1 : n_gen0
    for i = 1 : n_gen
        plot(f_impz{1}(1,:),amp_impz{idx}(i,:))
        hold on
    end
end
xlabel('frequency [Hz]')
ylabel('amplitude')
xlim([-0.1,3])
title('Individual spectrum')
hold off


%% Plot results

plt_switch = 1; % if 1 plot

if plt_switch == 1
    % rotor speed deviation (omega)
    figure
    for idx = 1 : n_gen0
        subplot(4,4,idx)
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
    
    % rotor angle
    figure
    for idx = 1 : n_gen0
        subplot(4,4,idx)
        for i = 1 : n_gen
            plot(t_impz{idx},ang_impz{idx}(i,:)-ang_impz{idx}(i,1),'DisplayName',strcat('Gen',num2str(i)));
            hold on
        end
        xlabel('time [s]')
        ylabel('speed deviation')
        title(strcat('Gen',num2str(idx)))
        if idx == 1
            legend('Location','best','NumColumns',2)
        end
        sgtitle('Rotor angle impz with different input loc')
        hold off
    end

else
    disp('No plot of data.')
end


% % angle
% figure
% for i = 1 : n_gen
%     plot(t,mac_ang(i,:),'DisplayName',strcat('Gen',num2str(i)));
%     hold on
% end
% xlabel('time [s]')
% ylabel('Angle [rad]')
% title('Rotor angle')
% legend('Location','best','NumColumns',2)
% hold off