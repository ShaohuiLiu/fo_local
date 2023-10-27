clear all
clear global
close % close graphics windows


%% Run TD with periodic input at target generator

td_switch = 0; % if 1 run TD, o.w. 0 just load previous results

ambient_input = 0; % if 1 add also ambient input 

dt_const = 0.005;  % uniform time step, o.w. use default in .m sys file.
T_tot = 150; % total time, if 0 then use default 15s

n_gen0 = 16; % simulate 1 fault location
fault_loc = 1:n_gen0;
% FO_input = [0.142857 0 0.005];% frequency, phase, amplitude local
% FO_input = [0.43 0 0.005];% frequency, phase, amplitude inter-area
% FO_input = [0.714 0 0.005];% frequency, phase, amplitude inter-area
% FO_input = [0.5 0 0.005];% frequency, phase, amplitude

FO_input = [0.57 0 0.005];% frequency, phase, amplitude no mode (0.57,0.93)

% dfile = 'data16m_ambient.m'; % pss with damping
dfile = 'data16m_damp_nopss1.m'; % 2 dampings, no pss, 03/14/23
pathname = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/';

simu_flag = zeros(n_gen0,1);
if td_switch == 1
    for i = 1 : n_gen0
        mac_fault = fault_loc(i);% fault location    
        % Run TD simulation
        simu_flag(i) = s_simu64_fo(pathname,dfile,mac_fault,FO_input,dt_const,T_tot,ambient_input);
    end
else
    simu_flag = ones(n_gen0,1);
end


%% Load simulation results
freq_impz = cell(n_gen0,1);
ang_impz = cell(n_gen0,1);
t_impz = cell(n_gen0,1);
n_gen = zeros(n_gen0,1); % number of generators
n_t = zeros(n_gen0,1); % number of samples

result_path = '/Users/shaohui/Documents/MATLAB/PST/68_bus_tests/simulation_result/';
for idx = 1 : n_gen0
    mac_fault = fault_loc(idx);
    if simu_flag(idx) == 1
        filename = strcat(result_path,'fo_td_fault',num2str(mac_fault),'.mat');
        load(filename)
        ang_impz{idx} = mac_ang;
        t_impz{idx} = t;
        n_gen(idx) = size(mac_ang,1); % n_gen in simulation results, overide
        disp(n_gen(idx))
        n_t(idx) = size(mac_ang,2);
        freq_impz{idx} = mac_spd - ones(n_gen(idx),n_t(idx)); % speed deviation
    else
        freq_impz{idx} = 'Nan'; % failed simulation
        ang_impz{idx} = 'Nan';
        t_impz{idx} = 'Nan';
        disp(strcat('Case ',num2str(idx),' failed!'))
    end
end

n_gen = ones(n_gen0,1) .* n_gen0; % n_gen in simulation results, overide

%% FFT analysis

% set t_start > 20 to remove unstable initial points
% t_start0 = find(t>20,1);
t_start0 = 1;

f_impz = cell(n_gen0,1);
amp_impz = cell(n_gen0,1);

% fft_band = [0.1,1];
for idx = 1 : n_gen0
%     [f_temp,p_temp] = spectrum_analysis(freq_impz{idx}(1,:),t_impz{idx},[]);
    [f_temp,p_temp] = spectrum_analysis(freq_impz{idx}(1,t_start0:end),t_impz{idx}(t_start0:end),[]);
    f_impz{idx} = zeros(n_gen(idx),length(f_temp));
    amp_impz{idx} = zeros(n_gen(idx),length(f_temp));
    for i = 1 : n_gen(idx)
        [f_impz{idx}(i,:),amp_impz{idx}(i,:)] = spectrum_analysis(freq_impz{idx}(i,t_start0:end),t_impz{idx}(t_start0:end),[]);
    end
end

%% plot spectrum
plt_switch = 1; % if 1 plot

if plt_switch == 1
    % fft of rotor speed deviation (omega)
    figure
    for idx = 1 : n_gen0
        subplot(ceil(sqrt(n_gen0)),floor(sqrt(n_gen0)),idx)
        for i = 1 : n_gen
            if i == fault_loc(idx)
                plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'LineWidth',2,'DisplayName',strcat('Gen',num2str(i),'(I)'));
                hold on
            else
                if max(amp_impz{idx}(i,:)) >= max(amp_impz{idx}(fault_loc(idx),:)) % oscillation >= source
                    plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'DisplayName',strcat('Gen',num2str(i),'*'));
                    hold on
                else
                    plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'LineStyle','-.','DisplayName',strcat('Gen',num2str(i)));
                    hold on
                end
            end
        end
        xlabel('frequency [Hz]')
        ylabel('amplitude')
        xlim([-0.5,2])
        title(strcat('Gen',num2str(idx)))
        if idx == 1
            legend('Location','best','NumColumns',2)
        end
        sgtitle('FFT of \omega FO responses')
        hold off
    end

else
    disp('No plot of spectrum.')
end

%% plot individual spectrum
plt_switch = 1; % if 1 plot

if plt_switch == 1
    % fft of rotor speed deviation (omega)
    figure
    for idx = 1 : n_gen0
        for i = 1 : n_gen
            subplot(4,4,i)
            if i == fault_loc(idx)
                plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'LineWidth',2,'DisplayName',strcat('Gen',num2str(i),'(I)'));
%                 hold on
                xlabel('frequency [Hz]')
                ylabel('amplitude')
                xlim([-0.5,4])
                title(strcat('Gen',num2str(i),'(',num2str(max(amp_impz{idx}(i,:))),')'))
%                 if idx == 1
%                 legend('Location','best','NumColumns',1)
%                 end
            else
                if max(amp_impz{idx}(i,:)) >= max(amp_impz{idx}(fault_loc(idx),:)) % oscillation >= source
                    plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'DisplayName',strcat('Gen',num2str(i),'*'));
                    hold on
                    plot(f_impz{idx}(fault_loc(idx),:),amp_impz{idx}(fault_loc(idx),:),'LineStyle','-.',...
                        'LineWidth',1,'DisplayName',strcat('Gen',num2str(fault_loc(idx)),'(I)'));
                else
                    plot(f_impz{idx}(i,:),amp_impz{idx}(i,:),'DisplayName',strcat('Gen',num2str(i)));
                    hold on
                    plot(f_impz{idx}(fault_loc(idx),:),amp_impz{idx}(fault_loc(idx),:),'LineStyle','-.',...
                        'LineWidth',1,'DisplayName',strcat('Gen',num2str(fault_loc(idx)),'(I)'));
                end
                xlabel('frequency [Hz]')
                ylabel('amplitude')
                xlim([-0.5,4])
                title(strcat('Gen',num2str(i),'(',num2str(max(amp_impz{idx}(i,:))),')'))
%                 if idx == 1
%                 legend('Location','best','NumColumns',1)
%                 end
            end
        end
        sgtitle('FFT of \omega FO responses')
        hold off
    end

else
    disp('No plot of spectrum.')
end

%%
idx = 16;
% figure
% 
% for i = 1 : 16
%     plot(t_impz{idx},freq_impz{idx}(i,:),'DisplayName',strcat('Gen',num2str(i)));
%     hold on
% end
% xlabel('time [s]')
% ylabel('speed deviation')
% title(strcat('Gen',num2str(idx)))
% if idx == 1
%     legend('Location','best','NumColumns',2)
% end
% hold off

figure
for i = 1 : 16
    if i == fault_loc(idx)
        plot(f_impz{idx}(i,2:end),amp_impz{idx}(i,2:end),'LineWidth',2,'DisplayName',strcat('Gen',num2str(i),'(I)'));
        hold on
    else
        if max(amp_impz{idx}(i,:)) >= max(amp_impz{idx}(fault_loc(idx),:)) % oscillation >= source
            plot(f_impz{idx}(i,2:end),amp_impz{idx}(i,2:end),'DisplayName',strcat('Gen',num2str(i),'*'));
            hold on
        else
            plot(f_impz{idx}(i,2:end),amp_impz{idx}(i,2:end),'LineStyle','-.','DisplayName',strcat('Gen',num2str(i)));
            hold on
        end
    end
end
xlabel('frequency [Hz]')
ylabel('amplitude')
xlim([-0.5,2])
title(strcat('Gen',num2str(idx)))


%% Plot results

plt_switch = 1; % if 1 plot

if plt_switch == 1
    % rotor speed deviation (omega)
    figure
    for idx = 1 : n_gen0
        subplot(ceil(sqrt(n_gen0)),floor(sqrt(n_gen0)),idx)
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
        subplot(ceil(sqrt(n_gen0)),floor(sqrt(n_gen0)),idx)
        for i = 1 : n_gen
            plot(t_impz{idx},ang_impz{idx}(i,:),'DisplayName',strcat('Gen',num2str(i)));
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

    % input signal
    noise0 = FO_input(3) .* cos((2.*pi.*FO_input(1)).*t_impz{idx} + FO_input(2));
    figure
    plot(t_impz{idx},noise0)
    xlabel('time [s]')
    ylabel('amplitude')
    title(strcat('Input at Gen',num2str(fault_loc(1))))

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