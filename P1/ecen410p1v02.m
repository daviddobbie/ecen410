%% ECEN 410 Project 1 - Distributed MIMO
% David Dobbie - 300340161




clear
clc

set(0,'defaultTextInterpreter','latex');



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





power_cap_results_no_dist = no_dist_MIMO(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R);

power_cap_results_dist_no_intercell = dist_MIMO_no_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R);



figure(2)
clf
hold on
p3 = plot(SNR_axis_dB,power_cap_results_no_dist);
p3.LineWidth = 2;
p3.LineStyle = '--'
p4 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell);
p4.LineWidth = 2;
p4.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist_{uplink}', 'Has dist_{uplink}')
lgnd.Location = 'NorthWest';





% This simulation holds pathloss constant and deals with pathloss, 
% shadowing, and Rayleigh fading. This holds as the BS antennas are 
% all in one location. This does not model inter-cell interference.
function power_cap_results_zf = no_dist_MIMO(uplink_pow_axis, K_users, ...
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

                L = db2pow(sigma_sfdB * randn(K_users,1));

                beta = (magnitude.^(-pathloss_m)) .* L;
                %beta = ones(K_users,1);

                H_rayleigh = (  sqrt(1/(2))* ...
                    (   randn(BS_ant,K_users) + 1j*randn(BS_ant,K_users) )    );

                R_sum_zf = 0;    

                % for each user in rnd uplink chan
                for user_per_chan = 1:K_users

                    h_k = H_rayleigh(:,user_per_chan);
                    beta_k = beta(user_per_chan);


                    h_not_k = H_rayleigh;
                    h_not_k(:,user_per_chan) = [];

                    beta_not_k = beta;
                    beta_not_k(user_per_chan) = [];

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
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R)

    user_test_count = 1000;
    channel_rnd_count = 100;

    power_cap_results_zf = zeros(length(uplink_pow_axis),1);

    %ant_pos = zeros(BS_ant,1) + 1j*zeros(BS_ant,1); 
    ant_bearing = linspace(0,2*pi,(BS_ant+1));
    ant_bearing(BS_ant+1) = []; %so we don't have two antenna together
    
    ant_magnitude_dist = (2/3)*R;

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

                H_rayleigh = (  sqrt(1/(2))* ...
                    (   randn(BS_ant,K_users) + 1j*randn(BS_ant,K_users) )    );

                R_sum_zf = 0;    

                % for each user in rnd uplink chan
                for user_per_chan = 1:K_users
                    
                    dist_from_ant = abs(pos(user_per_chan) - ant_pos);
                    
                    % assume equal power allocation to BS ant.
                    
                    % shadowing assumed to be independent due to dist
                    % between antennae - we get gain from this
                    L = db2pow(sigma_sfdB * randn(BS_ant,1));   
                    
                    
                    beta_k = (dist_from_ant'.^(-pathloss_m)) .* L;

                                        
                    
                    
                    


                    %---------- ZF compute SINR

                    known_interfere_term_zf = inv(H_rayleigh'*H_rayleigh);
                    
                    % normalised power dist
                    SINR_zf_kth = sum(uplink_power * beta_k) /...
                        (BS_ant * n_sigma^2 * known_interfere_term_zf(user_per_chan,user_per_chan));


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