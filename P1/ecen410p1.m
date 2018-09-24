%% ECEN 410 Project 1 - Distributed MIMO
% David Dobbie - 300340161




clear
clc

set(0,'defaultTextInterpreter','latex');



R = 1000;
pathloss_m = 3;
sigma_sfdB =1; %the fading
macro_userDens = 10;






base_station_location = [0 , 0i];


% generate the individual users








% uplink
BS_ant = 6;


n_sigma = 1;
K_users = 4

%uplink_pow_axis = linspace(db2pow(-20),db2pow(20),10);


uplink_pow_axis = logspace(-2,2,10);

user_test_count = 100;
channel_rnd_count = 100;




% for each different SNR
power_cap_results_mrc = zeros(length(uplink_pow_axis),1);
power_cap_results_mmse = zeros(length(uplink_pow_axis),1);
power_cap_results_zf = zeros(length(uplink_pow_axis),1);

for uplink_power_indx = 1:length(uplink_pow_axis)
    uplink_power = uplink_pow_axis(uplink_power_indx);
    
    
    
    
    R_user_results_mrc = zeros(user_test_count,1);
    R_user_results_mmse = zeros(user_test_count,1);
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
        
        R_results_mrc = zeros(channel_rnd_count,1);
        R_results_mmse = zeros(channel_rnd_count,1);  
        R_results_zf = zeros(channel_rnd_count,1);  
        
        
        
        % for each rnd channel
        for channel_test_num = 1:channel_rnd_count
            
            beta = (magnitude.^(-pathloss_m));
            beta = ones(K_users,1);
            
            H_rayleigh = (  sqrt(1/(2))* ...
                (   randn(BS_ant,K_users) + 1j*randn(BS_ant,K_users) )    );
            
            R_sum_mrc = 0;
            R_sum_mmse = 0;
            R_sum_zf = 0;          
            
            
            
            % for each user in rnd uplink chan
            for user_per_chan = 1:K_users
                
                h_k = H_rayleigh(:,user_per_chan);
                beta_k = beta(user_per_chan);
                
                
                h_not_k = H_rayleigh;
                h_not_k(:,user_per_chan) = [];
                
                beta_not_k = beta;
                beta_not_k(user_per_chan) = [];
                
                %---------- MRC compute SINR, R_sum
                top_term_mrc = uplink_power * beta_k * norm(h_k)^4;
                interfere_term_mrc = 0;
                for int_user_indx = 1:K_users-1
                    interfere_term_mrc =  interfere_term_mrc + ...
                        beta_not_k(int_user_indx)*abs(h_not_k(:,int_user_indx)'*h_k)^2; 
                end
                bot_term_mrc = n_sigma^2 * norm(h_k)^2 + uplink_power * interfere_term_mrc;               
                SINR_mrc_kth = top_term_mrc / bot_term_mrc;
                
                R_k_mrc = log2(1 + SINR_mrc_kth);
                R_sum_mrc = R_sum_mrc + R_k_mrc;
                
                
               
                
                %---------- MMSE compute SINR
                
                interfere_term_mmse = 0;
                for int_user_indx = 1:K_users-1
                    interfere_term_mmse =  interfere_term_mmse + ...
                        (h_not_k(:,int_user_indx)*h_not_k(:,int_user_indx)'); 
                end               
                
                bot_term_mmse = uplink_power * beta_k * ...
                    interfere_term_mmse + eye(BS_ant);
                SINR_mmse_kth = uplink_power * beta_k * h_k' *pinv(bot_term_mmse) * h_k;
                R_k_mmse = log2(1 + SINR_mmse_kth);
                R_sum_mmse = R_sum_mmse + R_k_mmse;                
                
                
                %---------- ZF compute SINR
                
                interfere_term_zf = inv(H_rayleigh'*H_rayleigh);
                SINR_zf_kth = uplink_power * beta_k /...
                    (n_sigma^2 * interfere_term_zf(user_per_chan,user_per_chan));
                
                
                R_k_zf = log2(1 + SINR_zf_kth);
                R_sum_zf = R_sum_zf + R_k_zf;                    
                
                

                
                
            end
            R_results_mrc(channel_test_num) = R_sum_mrc;
            R_results_mmse(channel_test_num) = R_sum_mmse;           
            R_results_zf(channel_test_num) = R_sum_zf;             
        end
        
        R_sum_mean_mrc = mean(R_results_mrc);
        R_user_results_mrc(user_test_num) = R_sum_mean_mrc;     
        
        R_sum_mean_mmse = mean(R_results_mmse);
        R_user_results_mmse(user_test_num) = R_sum_mean_mmse;      
        
        
        R_sum_mean_zf = mean(R_results_zf);
        R_user_results_zf(user_test_num) = R_sum_mean_zf;         
        
    end
    
    power_cap_results_mrc(uplink_power_indx) = mean(R_user_results_mrc);
    power_cap_results_mmse(uplink_power_indx) = mean(R_user_results_mmse);    
    power_cap_results_zf(uplink_power_indx) = mean(R_user_results_zf);    
    
end





figure(2)
clf
hold on
p1 = plot(pow2db(uplink_pow_axis),power_cap_results_mrc);
p1.LineWidth = 2;
p1.LineStyle = '-.'
p2 = plot(pow2db(uplink_pow_axis),power_cap_results_mmse);
p2.LineWidth = 2;
p3 = plot(pow2db(uplink_pow_axis),power_cap_results_zf);
p3.LineWidth = 2;
p3.LineStyle = '--'
hold off
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
ylim([0 35])
lgnd = legend('MRC_{uplink}','MMSE_{uplink}','ZF_{uplink}')
lgnd.Location = 'NorthWest';

