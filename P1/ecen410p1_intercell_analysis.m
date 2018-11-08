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


generate_interfering_users(20, ...
     1e3, 6);

sprintf('With Intercell Interference')




power_cap_results_dist_with_intercell_one_third = dist_MIMO_with_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 1/3);
sprintf('Antennas at R * 1/3 distribution simulation done')

power_cap_results_dist_with_intercell_mid = dist_MIMO_with_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 2/3);
sprintf('Antennas at R * 2/3 distribution simulation done')

power_cap_results_dist_with_intercell_on_edge = dist_MIMO_with_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 1);
sprintf('Antennas at R distribution simulation done')

power_cap_results_no_dist_with_intercell = no_dist_MIMO_v2_with_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R);
sprintf('Antennas not distributed simulation done')

sprintf('No Intercell Interference')

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
p4.LineStyle = '--'
p5 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_mid);
p5.LineWidth = 2;
p5.LineStyle = '--'
p6 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_on_edge);
p6.LineWidth = 2;
p6.LineStyle = '--'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist uplink', '$\frac{1}{3}R$  away dist uplink', ...
     '$\frac{2}{3}R$  away dist uplink', '$R$ away dist uplink')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])


figure(99)
clf
hold on
p3 = plot(SNR_axis_dB,power_cap_results_no_dist);
p3.LineWidth = 2;
p3.LineStyle = '--'
p4 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_one_third);
p4.LineWidth = 2;
p4.LineStyle = '--'
p5 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_mid);
p5.LineWidth = 2;
p5.LineStyle = '--'
colormap default

p3 = plot(SNR_axis_dB,power_cap_results_no_dist_with_intercell, 'Color', [0, 0.4470, 0.7410]);
p3.LineWidth = 2;
p3.LineStyle = '-.'
p4 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_one_third, 'Color', [0.8500, 0.3250, 0.0980]);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_mid, 'Color', 	[0.9290, 0.6940, 0.1250]);
p5.LineWidth = 2;
p5.LineStyle = '-.'

hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist. no intercell', '$\frac{1}{3}R$ dist no intercell', ...
     '$\frac{2}{3}R$ dist no intercell', ...
     'No dist. with intercell', '$\frac{1}{3}R$ dist with intercell', ...
     '$\frac{2}{3}R$ dist with intercell')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])




figure(3)
clf
hold on
p3 = plot(SNR_axis_dB,power_cap_results_no_dist_with_intercell);
p3.LineWidth = 2;
p3.LineStyle = '--'
p4 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_one_third);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_mid);
p5.LineWidth = 2;
p5.LineStyle = '-.'
p6 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_on_edge);
p6.LineWidth = 2;
p6.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist intercell', '$\frac{1}{3}R$  away dist intercell', ...
     '$\frac{2}{3}R$  away dist intercell', '$R$ away dist intercell')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])

figure(4)
clf
hold on
p4 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_one_third);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_one_third);
p5.LineWidth = 2;
p5.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('$\frac{1}{3}R$ uplink no intercell', ...
     '$\frac{1}{3}R$ uplink with intercell')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])

figure(5)
clf
hold on
p4 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_mid);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_mid);
p5.LineWidth = 2;
p5.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('$\frac{2}{3}R$ dist uplink no intercell', ...
     '$\frac{2}{3}R$ dist uplink with intercell')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])

figure(6)
clf
hold on
p4 = plot(SNR_axis_dB,power_cap_results_dist_no_intercell_on_edge);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_dist_with_intercell_on_edge);
p5.LineWidth = 2;
p5.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('$R$ dist uplink no intercell', ...
     '$R$ dist uplink with intercell')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])

figure(7)
clf
hold on
p4 = plot(SNR_axis_dB,power_cap_results_no_dist);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_no_dist_with_intercell);
p5.LineWidth = 2;
p5.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist uplink no intercell', ...
     'No dist uplink with intercell')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])

figure(7)
clf
hold on
p4 = plot(SNR_axis_dB,power_cap_results_no_dist);
p4.LineWidth = 2;
p4.LineStyle = '-.'
p5 = plot(SNR_axis_dB,power_cap_results_no_dist_with_intercell);
p5.LineWidth = 2;
p5.LineStyle = '-.'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist uplink no intercell', ...
     'No dist uplink with intercell')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([-50 20])




% This simulation holds pathloss constant and deals with pathloss, 
% shadowing, and Rayleigh fading. This holds as the BS antennas are 
% all in one location. This does not model inter-cell interference.
function power_cap_results_zf = no_dist_MIMO_v2(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R)

    user_test_count = 500;
    channel_rnd_count = 200;

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
                    SINR_zf_kth = uplink_power * BS_ant*  beta_k /...
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
% all in one location. This does not model inter-cell interference.
function power_cap_results_zf = no_dist_MIMO_v2_with_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R)

    user_test_count = 500;
    channel_rnd_count = 200;

    power_cap_results_zf = zeros(length(uplink_pow_axis),1);

    ant_pos = zeros(BS_ant,1);
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

            users_intercell_pos = generate_interfering_users(K_users, R, 6);

            % for each rnd channel
            for channel_test_num = 1:channel_rnd_count
                min_dim = min(K_users, BS_ant);
                H_rayleigh = (  sqrt(1/(2*min_dim))* ...
                    (   randn(BS_ant,K_users) + 1j*randn(BS_ant,K_users) )    );

                R_sum_zf = 0;    
                dist_from_ant = zeros(BS_ant, K_users);
                
                for user_per_chan = 1:K_users
                    dist_from_ant(:,user_per_chan) = abs(pos(user_per_chan) - ant_pos)';
                end

                % shadowing assumed to be independent due to dist
                % between antennae - we get gain from this
                L = db2pow(sigma_sfdB * randn(K_users,BS_ant));   
                %for where shadowing is the same for each antenna
                %L = db2pow(sigma_sfdB * randn(1,1)); 
                beta = ((dist_from_ant'.^(-pathloss_m)) .* L)';
                
                G = beta.*H_rayleigh; % channel with pl and shadow
                
                A = pinv(G); %the coder
                A = pinv(H_rayleigh); % only channel inversion
               
                
                
                % for each user in rnd uplink chan
                for user_per_chan = 1:K_users
                    
                    beta_k = beta(user_per_chan,1); % only one as we are dealing with no distribution

                    
                    
                    intercell_interfere = calculate_intercell_interference(uplink_power, ...
                    ant_pos, users_intercell_pos, sigma_sfdB, pathloss_m, A, user_per_chan);
                    %---------- ZF compute SINR

                    known_interfere_term_zf = inv(H_rayleigh'*H_rayleigh);
                    
                    SINR_zf_kth = uplink_power * BS_ant * beta_k /...
                        (n_sigma^2 * known_interfere_term_zf(user_per_chan,user_per_chan) + intercell_interfere);


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
function power_cap_results_zf = dist_MIMO_with_intercell(uplink_pow_axis, K_users, ...
    n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, dist_from_centre)

    user_test_count = 500;
    channel_rnd_count = 100;

    power_cap_results_zf = zeros(length(uplink_pow_axis),1);

    %ant_pos = zeros(BS_ant,1) + 1j*zeros(BS_ant,1); 
    ant_bearing = linspace(0,2*pi,(BS_ant+1));
    ant_bearing(BS_ant+1) = []; %so we don't have two antenna together
    
    ant_magnitude_dist = dist_from_centre*R;
    ant_pos = ant_magnitude_dist.* exp(1i*ant_bearing);
   
    
    for uplink_power_indx \= 1:length(uplink_pow_axis)
        uplink_power = uplink_pow_axis(uplink_power_indx);

        R_user_results_zf = zeros(user_test_count,1);    


        % on average users rnd gen
        for user_test_num = 1:user_test_count      
            % generate num of users
            magnitude = sqrt(abs(rand(K_users,1)*R^2));
            bearing = 2*pi*(rand(K_users,1));
            pos = magnitude .* exp(1i*bearing); 
            
            
            
            users_intercell_pos = generate_interfering_users(K_users, R, 6);
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
                
                dist_from_ant = zeros(BS_ant, K_users);
                
                for user_per_chan = 1:K_users
                    dist_from_ant(:,user_per_chan) = abs(pos(user_per_chan) - ant_pos)';
                end

                % shadowing assumed to be independent due to dist
                % between antennae - we get gain from this
                L = db2pow(sigma_sfdB * randn(K_users,BS_ant));   
                %for where shadowing is the same for each antenna
                %L = db2pow(sigma_sfdB * randn(1,1)); 
                beta = ((dist_from_ant'.^(-pathloss_m)) .* L)';
                
                G = beta.*H_rayleigh; % channel with pl and shadow
                
                A = pinv(G); %the coder
                A = pinv(H_rayleigh);
                % for each user in rnd uplink chan
                for user_per_chan = 1:K_users
                    
                    beta_k = beta(user_per_chan,:);

                    %---------- ZF compute SINR
                    known_interfere_term_zf = inv(H_rayleigh'*H_rayleigh);
                    
                    
                    intercell_interfere = calculate_intercell_interference(uplink_power, ...
                    ant_pos, users_intercell_pos, sigma_sfdB, pathloss_m, A, user_per_chan);
                    
                    % normalised power dist (assume equal power allocation)
                    SINR_zf_kth = uplink_power * 	sum(beta_k) /...
                        ( n_sigma^2 * known_interfere_term_zf(user_per_chan,user_per_chan) + intercell_interfere); %k,k element of known interfere
                    
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

% generates the interfering user positions relative to the centre of the
% primary cell we are tsting. As this is uplink, we only require the
% position of the interfereing users (assuming that all neighbouring cell
% have a synchronised uplink session).

function user_pos = generate_interfering_users(K_users, ...
     R, num_neighbour_cells)
 
    cell_bearing = linspace(0,2*pi,(num_neighbour_cells+1));
    cell_bearing(num_neighbour_cells+1) = []; %so we don't have two antenna together
    
    cell_magnitude_dist = 2*R;
    cell_pos = cell_magnitude_dist.* exp(1i*cell_bearing);
    
    magnitude = sqrt(abs(rand(K_users*num_neighbour_cells,1)*R^2));
    bearing = 2*pi*(rand(K_users*num_neighbour_cells,1));
    user_pos = magnitude .* exp(1i*bearing); 
    
    for cell_indx = 1:num_neighbour_cells
        top = K_users*(cell_indx-1) + 1;
        bot = K_users*(cell_indx);
        user_pos(top:bot) = user_pos(top:bot) + cell_pos(cell_indx) ;
    end
%{
    figure(66)
    clf
    
    figure(54)
    clf
    hold on
    scatter(real(user_pos),imag(user_pos),'b.')
    scatter(0,0,'r+')
    hold off
    xlim([-3*R,3*R])
    ylim([-3*R,3*R])
    grid on  
%}
end

function power = calculate_intercell_interference(uplink_power, ...
    antenna_pos, user_pos, sigma_sfdB, pathloss_m, A, user_num)

    BS_ant = length(antenna_pos);
    
    ant_inter = zeros(BS_ant,1);
    
    for antenna_idx = 1:BS_ant
        antenna = antenna_pos(antenna_idx);

        dist_from_ant = abs(user_pos - antenna);
        
        % slow fading
        L = db2pow(sigma_sfdB * randn(length(user_pos),1)); 
        beta = (dist_from_ant'.^(-pathloss_m))' .* L;        
        h_rayleigh = (sqrt(1/(2))*  ( randn(1) + 1j*randn(1) ));
        
        ant_inter(antenna_idx) = A(user_num,antenna_idx)*uplink_power*sum((beta*h_rayleigh));
        
    end
    
    power = sum(abs(ant_inter));
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
                    SINR_zf_kth = uplink_power * sum(beta_k) /...
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