%% ECEN 410 Project 1 - Distributed MIMO
% David Dobbie - 300340161




clear
clc

set(groot,'defaultLineLineWidth',2)
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)


R = 1000;
pathloss_m = 4.6;
sigma_sfdB =10; %the shadowing

% uplink
BS_ant = 6;

boltzmann = 1.3807e-23;
ambient_temp = 290; %20 celsius
noise_bandwidth = 5e6; %5 Mhz of WCDMA channel
n_sigma = boltzmann * ambient_temp * noise_bandwidth;
K_users = 4


base_station_location = [0 , 0i];
%uplink_pow_axis = linspace(db2pow(-20),db2pow(20),10);



SNR_axis_dB = pow2db(logspace(-5,2,10));

uplink_pow_axis =db2pow(SNR_axis_dB) * n_sigma;





power_cap_results_no_dist = no_dist_MIMO_v2(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R);

sprintf('Antennas not distributed simulation done')

power_cap_results_dist_no_intercell_one_third = dist_MIMO_no_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 1/3);
sprintf('Antennas at R * 1/3 distribution simulation done')

power_cap_results_dist_no_intercell_mid = dist_MIMO_no_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 2/3);
sprintf('Antennas at R * 2/3 distribution simulation done')

power_cap_results_dist_no_intercell_on_edge = dist_MIMO_no_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 1);
sprintf('Antennas at R distribution simulation done')


figure(2)
clf
hold on
p3 = plot(SNR_axis_dB,power_cap_results_no_dist);
p3.LineWidth = 2;
p3.LineStyle = '--'
p4 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_one_third);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_mid);
p5.LineWidth = 2;
p5.LineStyle = '-.'
p6 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_on_edge);
p6.LineWidth = 2;
p6.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist uplink', '$\frac{1}{3}R$  away dist uplink', ...
     '$\frac{2}{3}R$  away dist uplink', '$R$ away dist uplink')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';


% This simulation holds pathloss constant and deals with pathloss, 
% shadowing, and Rayleigh fading. This holds as the BS antennas are 
% all in one location. This does not model inter-cell interference.
function power_cap_results_zf = no_dist_MIMO_v2(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R)

    user_test_count = 50;
    channel_rnd_count = 100;

    power_cap_results_zf = zeros(length(uplink_pow_axis),1);

    for uplink_power_indx = 1:length(uplink_pow_axis)
        uplink_power = uplink_pow_axis(uplink_power_indx);

        R_user_results_zf = zeros(user_test_count,1);    


        % on average users rnd gen
        for user_test_num = 1:user_test_count      
            % generate num of users
            magnitude = sqrt(abs(rand(K_users,1)*R^2));
            bearing = 2*pi*(rand(K_users,1));
            pos = magnitude .* exp(1i*bearing);    


            %----------------
            %{
            figure(1)
            clf
            hold on
            scatter(real(pos),imag(pos),'b.')
            scatter(base_station_location(1),base_station_location(2),'r+')
            hold off
            xlim([-R,R])
            ylim([-R,R])
            grid on  
         %}
            %-------------------

            R_results_zf = zeros(channel_rnd_count,1);  



            % for each rnd channel
            for channel_test_num = 1:channel_rnd_count
                min_dim = min(K_users, BS_ant);
                H_rayleigh = (  sqrt(1/(2*min_dim))* ...
                    (   randn(BS_ant,K_users) + 1j*randn(BS_ant,K_users) )    );

                R_sum_zf = 0;    

                % for each user in rnd uplink chan
                for user_per_chan = 1:K_users
                    
                    dist_from_ant = abs(magnitude(user_per_chan));
                    
                    
                    % shadowing assumed to be independent due to dist
                    % between antennae - we get gain from this
                    L = db2pow(sigma_sfdB * randn(1,1));   
                    
                    
                    %L = db2pow(sigma_sfdB * randn(1,1));
                    
                    
                    beta_k = (dist_from_ant'.^(-pathloss_m)) .* L;
                    
                    %h_k = H_rayleigh(:,user_per_chan);
                    %beta_k = beta(user_per_chan);


                    %h_not_k = H_rayleigh;
                    %h_not_k(:,user_per_chan) = [];

                    %beta_not_k = beta;
                    %beta_not_k(user_per_chan) = [];

                    %---------- ZF compute SINR

                    known_interfere_term_zf = inv(H_rayleigh'*H_rayleigh);
                    SINR_zf_kth = uplink_power * beta_k /...
                        (n_sigma^2 * known_interfere_term_zf(user_per_chan,user_per_chan));


                    R_k_zf = log2(1 + SINR_zf_kth);
                    R_sum_zf = R_sum_zf + R_k_zf;                    

                end         
                R_results_zf(channel_test_num) = R_sum_zf;             
            end

            R_sum_mean_zf = mean(R_results_zf);
            R_user_results_zf(user_test_num) = R_sum_mean_zf;         

        end

        power_cap_results_zf(uplink_power_indx) = mean(R_user_results_zf);    

    end

end

















% This simulation holds pathloss constant and deals with pathloss, 
% shadowing, and Rayleigh fading. This holds as the BS antennas are 
% in different locations. This does not model inter-cell interference.
function power_cap_results_zf = dist_MIMO_no_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, dist_from_centre)

    user_test_count = 500;
    channel_rnd_count = 100;

    power_cap_results_zf = zeros(length(uplink_pow_axis),1);

    %ant_pos = zeros(BS_ant,1) + 1j*zeros(BS_ant,1); 
    ant_bearing = linspace(0,2*pi,(BS_ant+1));
    ant_bearing(BS_ant+1) = []; %so we don't have two antenna together
    
    ant_magnitude_dist = dist_from_centre*R;
    ant_pos = ant_magnitude_dist.* exp(1i*ant_bearing);
   
    
    for uplink_power_indx = 1:length(uplink_pow_axis)
        uplink_power = uplink_pow_axis(uplink_power_indx);

        R_user_results_zf = zeros(user_test_count,1);    


        % on average users rnd gen
        for user_test_num = 1:user_test_count      
            % generate num of users
            magnitude = sqrt(abs(rand(K_users,1)*R^2));
            bearing = 2*pi*(rand(K_users,1));
            pos = magnitude .* exp(1i*bearing);    
            %----------------
            %{
            figure(1)
            clf
            hold on
            scatter(real(pos),imag(pos),'b.')
            scatter(base_station_location(1),base_station_location(2),'r+')
            hold off
            xlim([-R,R])
            ylim([-R,R])
            grid on  
         %}
            %-------------------
            R_results_zf = zeros(channel_rnd_count,1);  
            % for each rnd channel
            for channel_test_num = 1:channel_rnd_count
               %beta = ones(K_users,1);
                min_dim = min(K_users, BS_ant);
                H_rayleigh = (  sqrt(1/(2*min_dim))* ...
                    (   randn(BS_ant,K_users) + 1j*randn(BS_ant,K_users) )    );
                R_sum_zf = 0;    
                % for each user in rnd uplink chan
                for user_per_chan = 1:K_users
                    
                    dist_from_ant = abs(pos(user_per_chan) - ant_pos);
                    
                    
                    % shadowing assumed to be independent due to dist
                    % between antennae - we get gain from this
                    L = db2pow(sigma_sfdB * randn(BS_ant,1));   
                    
                    %for where shadowing is the same for each antenna
                    %L = db2pow(sigma_sfdB * randn(1,1)); 
                    
                    
                    
                    beta_k = (dist_from_ant'.^(-pathloss_m)) .* L;

                    %---------- ZF compute SINR
                    known_interfere_term_zf = inv(H_rayleigh'*H_rayleigh);
                    
                    % normalised power dist (assume equal power allocation)
                    SINR_zf_kth = sum(uplink_power/BS_ant * beta_k) /...
                        ( n_sigma^2 * known_interfere_term_zf(user_per_chan,user_per_chan)); %k,k element of knwon interfere
                    
                    R_k_zf = log2(1 + SINR_zf_kth);
                    R_sum_zf = R_sum_zf + R_k_zf;                    

                end         
                R_results_zf(channel_test_num) = R_sum_zf;             
            end

            R_sum_mean_zf = mean(R_results_zf);
            R_user_results_zf(user_test_num) = R_sum_mean_zf;         

        end

        power_cap_results_zf(uplink_power_indx) = mean(R_user_results_zf);    

    end

end