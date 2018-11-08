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



base_station_location = [0 , 0i];
%uplink_pow_axis = linspace(db2pow(-20),db2pow(20),10);



SNR_axis_dB = pow2db(20);

uplink_pow_axis =db2pow(SNR_axis_dB) * n_sigma;

K_user_axis = 1:1:6;





sprintf('With Intercell Interference')

zero_third_sum_rate = zeros(length(K_user_axis),1);
one_third_sum_rate = zeros(length(K_user_axis),1);
two_third_sum_rate = zeros(length(K_user_axis),1);
three_third_sum_rate = zeros(length(K_user_axis),1);



for indx_users = 1:length(K_user_axis)
    K_users = K_user_axis(indx_users)
    
    one_third_sum_rate(indx_users) = dist_MIMO_with_intercell(uplink_pow_axis, K_users, ...
        n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 1/3);
    sprintf('Antennas at R * 1/3 distribution simulation done')

    two_third_sum_rate(indx_users) = dist_MIMO_with_intercell(uplink_pow_axis, K_users, ...
        n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 2/3);
    sprintf('Antennas at R * 2/3 distribution simulation done')

    three_third_sum_rate(indx_users) = dist_MIMO_with_intercell(uplink_pow_axis, K_users, ...
        n_sigma, BS_ant, sigma_sfdB, pathloss_m, R, 1);
    sprintf('Antennas at R distribution simulation done')

    zero_third_sum_rate(indx_users) = no_dist_MIMO_v2_with_intercell(uplink_pow_axis, K_users, ...
        n_sigma, BS_ant, sigma_sfdB, pathloss_m, R);
    sprintf('Antennas not distributed simulation done')

end


figure(1)
clf
hold on
p3 = plot(K_user_axis,zero_third_sum_rate);
p3.LineWidth = 2;
p3.LineStyle = '--'
p3.Marker='+'
p4 = plot(K_user_axis,one_third_sum_rate);
p4.LineWidth = 2;
p4.LineStyle = '--'
p4.Marker='o'
p5 = plot(K_user_axis,two_third_sum_rate);
p5.LineWidth = 2;
p5.LineStyle = '--'
p5.Marker='square'
p6 = plot(K_user_axis,three_third_sum_rate);
p6.LineWidth = 2;
p6.LineStyle = '--'
p6.Marker='diamond'
hold off
grid on
xlabel('Users per cell ($K$)')
ylabel('Sum Rate (bits/s/Hz)')
%ylim([0 35])
lgnd = legend('No dist uplink', '$\frac{1}{3}R$  away dist uplink', ...
     '$\frac{2}{3}R$  away dist uplink', '$R$ away dist uplink')
lgnd.Location = 'NorthWest';
lgnd.Interpreter = 'latex';
xlim([1 6])


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
   
    
    for uplink_power_indx = 1:length(uplink_pow_axis)
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




