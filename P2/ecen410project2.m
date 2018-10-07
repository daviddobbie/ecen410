%% ECEN 410 - Project 2 - Spatial Channel Models
% David Dobbie

% AIM: to comapre the predicted single-user MIMO capacity between classical
% Kronecker correlated fading and ray-based channel modelling.

% EXPERIMENTAL DETAILS:
% We will examine the different ergodic MIMO capacity for no channel state
% information at the transmitter (no CSIT). The channel model here deals
% with the scattering and fading of the model. We will know that the
% ray based model is accurately reproduced if wil la alarge number of
% clusters it converges towards the Kronecker model.

clc
clear
close all

SNR_dB = 10;
dim = 4;



rho = 0.90;
ant = 0:1:dim-1;

phase_rho = exp(1j*2*pi*rand(4,4));
ant_exponent = toeplitz(ant);
Rx_corr = db2pow(SNR_dB)*phase_rho.^ant_exponent;


H_trials = 1e4;
capacity_classical = zeros(1,H_trials);  

for cap_idx = 1:H_trials
    % create random unit var, zero mean, normally distirbuted channel
    % normalise across 4 dimensional Gaussian to unit variance
    H_classical = sqrt(1/(2*dim))*(randn(dim) + 1j*randn(dim));
    term_classical = eye(dim) + H_classical*Rx_corr*ctranspose(H_classical);
    % add absolute term to make the logarithm behave
    capacity_classical(cap_idx) = log2(abs(det(term_classical))); 
end

ergodic_cap_mean = mean(capacity_classical);

figure(1)
clf
hold on
p1 = cdfplot(capacity_classical);
p1.LineWidth = 2;

%xlim([1 10])
ylabel('CDF')
xlabel('Rate (bps/Hz)')
title('')
lgnd = legend('Classical Fading Model');
grid on
hold off


%% Cluster Based Modelling

H_trials = 1e4;
capacity_cluster = zeros(1,H_trials);  

rays_per_cluster = 5;
cluster_count = 30;
sigma_s = 1;
sigma_c = 1;
SNR_dB = 10;
dim = 4;
wavelength = 3e8/(6e6);
pathloss_cluster  = 4;
cluster_dist_from_rx = 10;


for cap_idx = 1:H_trials
    
    capacity_cluster(cap_idx) = cluster_based_model(rays_per_cluster, ...
    cluster_count, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
    pathloss_cluster, cluster_dist_from_rx); 
end


cluster_based_model(rays_per_cluster, ...
    cluster_count, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
    pathloss_cluster, cluster_dist_from_rx)
%{
SNR_dB = 10;
dim = 4;
Rx_corr = db2pow(SNR_dB)*eye(dim);

wavelength = 3e8/(51e6);




rays_per_cluster = 50;
cluster_count = 1;
sigma_s = 1;
sigma_c = 1;

antenna_dist= wavelength/2;



antenna_dist = 1*antenna_dist;

base_station_coords = [0;0]; % base station at 0,0


pathloss_cluster = 3;
cluster_dist_from_rx = 10;

beta_cluster = cluster_dist_from_rx^(-pathloss_cluster);

ant = 0:1:dim-1;

for cap_idx = 1:H_trials
    
    H = [];
    summed_cluster_terms = zeros(dim,dim);
    
    
    for cluster_indx = 1:cluster_count
        cluster_coeff_phase = 2*pi*rand(cluster_count,1);
        complex_cluster_coeff = sqrt(beta_cluster) * exp(1j* cluster_coeff_phase);
        central_cluster_angle = sigma_c*randn(1, cluster_count);    
        

        summed_ray_terms = zeros(dim,dim);
        
        
        for ray_indx = 1:rays_per_cluster

            offset_angle = rand_laplace(0,1,rays_per_cluster);

            phi_AOD = central_cluster_angle(cluster_indx) + offset_angle(ray_indx);
            phi_AOA = central_cluster_angle(cluster_indx) + offset_angle(ray_indx);

            a_tx_AOD = exp(antenna_dist*j*2*pi*ant* cos(phi_AOD))';
            
            a_rx_AOA = exp(antenna_dist*j*2*pi*ant* cos(phi_AOA))';           
            
            h_iid_ray = (1\sqrt(2))*(randn(1,1) + j*randn(1,1));
            %h_iid_ray = 1;
            
            ray_term = h_iid_ray * a_rx_AOA * a_tx_AOD';
            summed_ray_terms = sqrt(complex_cluster_coeff(cluster_indx) ...
            /rays_per_cluster) * ray_term + summed_ray_terms;
            

            
        end
        %summed_cluster_terms = summed_ray_terms./cluster_count + summed_cluster_terms;
        summed_cluster_terms = summed_ray_terms + summed_cluster_terms;
        
        
    end
    H_cluster = sqrt(db2pow(SNR_dB))*summed_cluster_terms;
       
    term = eye(dim) + H_cluster*Rx_corr*ctranspose(H_cluster);
    
    % add absolute term to make the logarithm behave    
    capacity_cluster(cap_idx) = log2(abs(det(term))); 
end
%}



ergodic_cap_mean = mean(capacity_cluster);


figure(1)
hold on
p1 = cdfplot(capacity_cluster)
p1.LineWidth = 2;

%xlim([1 10])
ylabel('CDF')
xlabel('Rate (bps/Hz)')
title('')
lgnd = legend('Classical Fading Model','Cluster Fading Model')
grid on
hold off


function y = rand_laplace(mu, sigma, m)    
    u = rand(m, 1)-0.5;
    b = sigma / sqrt(2);
    y = mu - b * sign(u).* log(1- 2* abs(u));
end



function capacity_cluster = cluster_based_model(rays_per_cluster, ...
    cluster_count, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
    pathloss_cluster, cluster_dist_from_rx)

    Rx_corr = db2pow(SNR_dB)*eye(dim);
    antenna_dist= wavelength/2;

    beta_cluster = cluster_dist_from_rx^(-pathloss_cluster);

    ant = 0:1:dim-1;

    H = [];
    summed_cluster_terms = zeros(dim,dim);


    for cluster_indx = 1:cluster_count
        cluster_coeff_phase = 2*pi*rand(cluster_count,1);
        complex_cluster_coeff = sqrt(beta_cluster) * exp(1j* cluster_coeff_phase);
        central_cluster_angle = sigma_c*randn(1, cluster_count);    


        summed_ray_terms = zeros(dim,dim);


        for ray_indx = 1:rays_per_cluster

            offset_angle = rand_laplace(0,sigma_s,rays_per_cluster);

            phi_AOD = central_cluster_angle(cluster_indx) + offset_angle(ray_indx);
            phi_AOA = central_cluster_angle(cluster_indx) + offset_angle(ray_indx);

            a_tx_AOD = exp(antenna_dist*j*2*pi*ant* cos(phi_AOD))';
            a_rx_AOA = exp(antenna_dist*j*2*pi*ant* cos(phi_AOA))';           

            h_iid_ray = (1\sqrt(2))*(randn(1,1) + j*randn(1,1));
            %h_iid_ray = 1;

            ray_term = h_iid_ray * a_rx_AOA * a_tx_AOD';
            summed_ray_terms = sqrt(complex_cluster_coeff(cluster_indx) ...
            /rays_per_cluster) * ray_term + summed_ray_terms;



        end
        
        % normalise for cluster count
        summed_cluster_terms = summed_ray_terms./cluster_count + summed_cluster_terms;
        %summed_cluster_terms = summed_ray_terms + summed_cluster_terms;


    end
    H_cluster = sqrt(db2pow(SNR_dB))*summed_cluster_terms;

    term = eye(dim) + H_cluster*Rx_corr*ctranspose(H_cluster);

    % add absolute term to make the logarithm behave    
    capacity_cluster = log2(abs(det(term))); 

end

























