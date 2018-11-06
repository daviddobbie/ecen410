%% ECEN 410 - Project 2 - Spatial Channel Models
% David Dobbie

% AIM: to compare the predicted single-user MIMO capacity between classical
% Kronecker correlated fading and ray-based channel modelling.

% EXPERIMENTAL DETAILS:
% We will examine the different ergodic MIMO capacity for no channel state
% information at the transmitter (no CSIT). The channel model here deals
% with the scattering and fading of the model. We will know that the
% ray based model is accurately reproduced if wil la alarge number of
% clusters it converges towards the Kronecker model.

clc
clear

set(groot,'defaultLineLineWidth',2)

set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)





%% -------- Analysis of Ergodic Capacity Between Models
dim = 4;
H_trials = 1e3;

capacity_classical_corr1 = zeros(1,H_trials);  
capacity_classical_independ = zeros(1,H_trials);  

SNR_axis = 0:1:20;
%SNR_axis = 10;
ergodic_cap_results_classical_corr1 = zeros(1,length(SNR_axis));
ergodic_cap_results_classical_independ = zeros(1,length(SNR_axis));

for SNR_indx = 1:length(SNR_axis)
    SNR_dB = SNR_axis(SNR_indx)
    

    
    for cap_idx = 1:H_trials    
        capacity_classical_corr1(cap_idx) = kronecker_based_model(SNR_dB, dim, 0.95);
        capacity_classical_independ(cap_idx) = kronecker_based_model(SNR_dB, dim, 0);
    end
    %ergodic_cap_mean = mean(capacity_classical);
    ergodic_cap_results_classical_corr1(SNR_indx) = mean(capacity_classical_corr1);
    ergodic_cap_results_classical_independ(SNR_indx) = mean(capacity_classical_independ);
    
    figure(1)
    clf
    hold on
    p1 = cdfplot(capacity_classical_corr1);
    p2 = cdfplot(capacity_classical_independ);

    %xlim([1 10])
    ylabel('CDF')
    xlabel('Rate (bps/Hz)')
    title('')
    grid on
    hold off
end

figure(2)
clf
hold on
plot(SNR_axis,ergodic_cap_results_classical_corr1);
plot(SNR_axis,ergodic_cap_results_classical_independ);
legend('Kronecker \rho = 1', 'Kronecker \rho = 0');

% Cluster Based Modelling



rays_per_cluster = 20;
cluster_count = 20;
% sigma_s = 2*pi*(5/360);
% sigma_c = 2*pi*(11/360);
sigma_s = 2*pi*(5/360);
sigma_c = 0;
SNR_dB = 10;
dim = 4;
wavelength = 3e8/(2.4e9);
pathloss_cluster  = 4.58;
cluster_dist_from_rx = 20;
cluster_shadowing = 3;

%SNR_axis = 0:1:20;
%SNR_axis = 10;

capacity_cluster = zeros(1,H_trials);  
ergodic_cap_results_cluster = zeros(1,length(SNR_axis));

for SNR_indx = 1:length(SNR_axis)
    tic
    SNR_dB = SNR_axis(SNR_indx)
    
    for cap_idx = 1:H_trials

        capacity_cluster(cap_idx) = cluster_based_model(rays_per_cluster, ...
        cluster_count, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 

    end


    ergodic_cap_results_cluster(SNR_indx) = mean(capacity_cluster);

    figure(1)
    hold on
    p1 = cdfplot(capacity_cluster);
    p1.LineWidth = 2;
    hold off
    %xlim([1 10])
    ylabel('CDF');
    xlabel('Rate (bps/Hz)');
    title('');
    lgnd = legend('Classical Fading Model','Cluster Fading Model');
    grid on
    toc
end

figure(2)
hold on
p3 = plot(SNR_axis,ergodic_cap_results_cluster);
p3.LineWidth = 2;
hold off
legend('Kronecker \rho = 0.95', 'Kronecker \rho = 0','Cluster Model');
ylabel('Ergodic Capcacity (bps/Hz)')
xlabel('SNR(dB)')
grid on


%% --- Analysis of CDF 10% and 1% cutoff

dim = 4;
H_trials = 1e3;

capacity_classical_corr1 = zeros(1,H_trials);  
capacity_classical_independ = zeros(1,H_trials);  

SNR_axis = 0:1:20;
%SNR_axis = 10;
cdf_A_cap_results_classical_corr1 = zeros(1,length(SNR_axis));
cdf_A_cap_results_classical_independ = zeros(1,length(SNR_axis));

cdf_B_cap_results_classical_corr1 = zeros(1,length(SNR_axis));
cdf_B_cap_results_classical_independ = zeros(1,length(SNR_axis));

for SNR_indx = 1:length(SNR_axis)
    SNR_dB = SNR_axis(SNR_indx)
    tic

    
    for cap_idx = 1:H_trials    
        capacity_classical_corr1(cap_idx) = kronecker_based_model(SNR_dB, dim, 0.95);
        capacity_classical_independ(cap_idx) = kronecker_based_model(SNR_dB, dim, 0);
    end
    
    [f ,x] = ecdf(capacity_classical_corr1);
    [c idx] = min(abs(f - 0.1));
    cdf_A_corr1 = x(idx);
    
    [c idx] = min(abs(f - 0.05));
    cdf_B_corr1 = x(idx);
    
    
    
    [f ,x] = ecdf(capacity_classical_independ);
    [c idx] = min(abs(f - 0.1)); 
    cdf_A_corr_indep = x(idx);    
    
    [c idx] = min(abs(f - 0.05));
    cdf_B_corr_indep = x(idx);   
    
    %ergodic_cap_mean = mean(capacity_classical);
    cdf_A_cap_results_classical_corr1(SNR_indx) = cdf_A_corr1;
    cdf_A_cap_results_classical_independ(SNR_indx) = cdf_A_corr_indep;
    
    cdf_B_cap_results_classical_corr1(SNR_indx) = cdf_B_corr1;
    cdf_B_cap_results_classical_independ(SNR_indx) = cdf_B_corr_indep;
    
    
    figure(1)
    clf
    hold on
    p1 = cdfplot(capacity_classical_corr1);
    p2 = cdfplot(capacity_classical_independ);

    %xlim([1 10])
    ylabel('CDF')
    xlabel('Rate (bps/Hz)')
    title('')
    grid on
    hold off
    toc
end

figure(3)
clf
hold on
plot(SNR_axis,cdf_A_cap_results_classical_corr1);
plot(SNR_axis,cdf_A_cap_results_classical_independ);
plot(SNR_axis,cdf_B_cap_results_classical_corr1);
plot(SNR_axis,cdf_B_cap_results_classical_independ);
legend('Kronecker \rho = 1, 10% Outage', 'Kronecker \rho = 0 10% Outage' ,...
    'Kronecker \rho = 1, 5% Outage', 'Kronecker \rho = 0 5% Outage');

% Cluster Based Modelling



rays_per_cluster = 20;
cluster_count = 20;
% sigma_s = 2*pi*(5/360);
% sigma_c = 2*pi*(11/360);
sigma_s = 2*pi*(5/360);
sigma_c = 0;
SNR_dB = 10;
dim = 4;
wavelength = 3e8/(2.4e9);
pathloss_cluster  = 4.58;
cluster_dist_from_rx = 20;
cluster_shadowing = 3;

%SNR_axis = 0:1:20;
%SNR_axis = 10;

capacity_cluster = zeros(1,H_trials);  

cdf_A_cap_results_cluster = zeros(1,length(SNR_axis));
cdf_B_cap_results_cluster = zeros(1,length(SNR_axis));


for SNR_indx = 1:length(SNR_axis)
    tic
    SNR_dB = SNR_axis(SNR_indx)
    
    for cap_idx = 1:H_trials

        capacity_cluster(cap_idx) = cluster_based_model(rays_per_cluster, ...
        cluster_count, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 

    end

    [f ,x] = ecdf(capacity_cluster);
    [c idx] = min(abs(f - 0.1)); 
    cdf_A_cluster = x(idx);    
    
    [c idx] = min(abs(f - 0.05));
    cdf_B_cluster = x(idx);   
    
    %ergodic_cap_mean = mean(capacity_classical);
    cdf_A_cap_results_cluster(SNR_indx) = cdf_A_cluster;
    cdf_B_cap_results_cluster(SNR_indx) = cdf_B_cluster;


    figure(1)
    hold on
    p1 = cdfplot(capacity_cluster);
    p1.LineWidth = 2;
    hold off
    %xlim([1 10])
    ylabel('CDF');
    xlabel('Rate (bps/Hz)');
    title('');
    lgnd = legend('Classical Fading Model','Cluster Fading Model');
    grid on
    toc
end

figure(3)
hold on
p3 = plot(SNR_axis,cdf_A_cap_results_cluster);
p3.LineWidth = 2;
p3 = plot(SNR_axis,cdf_B_cap_results_cluster);
p3.LineWidth = 2;
hold off
legend('Kronecker \rho = 1, 10% Outage', 'Kronecker \rho = 0, 10% Outage' ,...
    'Kronecker \rho = 1, 5% Outage', 'Kronecker \rho = 0, 5% Outage', ...
    'Cluster, 10% Outage', 'Cluster, 5% Outage');
ylabel('Outage Capacity (bps/Hz)')
xlabel('SNR(dB)')
grid on



%% --- Analysis of Different Cluster Types

%


rays_per_cluster = 20;
cluster_count = 20;
% sigma_s = 2*pi*(5/360);
% sigma_c = 2*pi*(11/360);
sigma_s = 2*pi*(5/360);
sigma_c = 0;
SNR_dB = 10;
dim = 4;
wavelength = 3e8/(2.4e9);
pathloss_cluster  = 4.58;
cluster_dist_from_rx = 20;
cluster_shadowing = 3;

%SNR_axis = 0:1:20;
%SNR_axis = 10;

capacity_cluster_test_A = zeros(1,H_trials);  
capacity_cluster_test_B = zeros(1,H_trials);  
capacity_cluster_test_C = zeros(1,H_trials);  
capacity_cluster_test_D = zeros(1,H_trials);  

capacity_cluster_test_E = zeros(1,H_trials);  
capacity_cluster_test_F = zeros(1,H_trials);  
capacity_cluster_test_G = zeros(1,H_trials);  
capacity_cluster_test_H = zeros(1,H_trials);  

capacity_cluster_test_I = zeros(1,H_trials);  



ergodic_cap_test_A = zeros(1,length(SNR_axis));
ergodic_cap_test_B = zeros(1,length(SNR_axis));
ergodic_cap_test_C = zeros(1,length(SNR_axis));
ergodic_cap_test_D = zeros(1,length(SNR_axis));

ergodic_cap_test_E = zeros(1,length(SNR_axis));
ergodic_cap_test_F = zeros(1,length(SNR_axis));
ergodic_cap_test_G = zeros(1,length(SNR_axis));
ergodic_cap_test_H = zeros(1,length(SNR_axis));

ergodic_cap_test_I = zeros(1,length(SNR_axis));  

for SNR_indx = 1:length(SNR_axis)
    tic
    SNR_dB = SNR_axis(SNR_indx)
    
    for cap_idx = 1:H_trials
        capacity_cluster_test_A(cap_idx) = cluster_based_model(20, ...
        20, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 
    
        capacity_cluster_test_B(cap_idx) = cluster_based_model(10, ...
        20, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 
    
        capacity_cluster_test_C(cap_idx) = cluster_based_model(5, ...
        20, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 
    
        capacity_cluster_test_D(cap_idx) = cluster_based_model(1, ...
        20, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 
    
        capacity_cluster_test_E(cap_idx) = capacity_cluster_test_A(cap_idx) ; % same var 
    
        capacity_cluster_test_F(cap_idx) = cluster_based_model(20, ...
        10, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 
    
        capacity_cluster_test_G(cap_idx) = cluster_based_model(20, ...
        5, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing); 
    
        capacity_cluster_test_H(cap_idx) = cluster_based_model(20, ...
        1, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing);  

        capacity_cluster_test_I(cap_idx) = cluster_based_model(1, ...
        1, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
        pathloss_cluster, cluster_dist_from_rx, cluster_shadowing);

    end
    ergodic_cap_test_A(SNR_indx) = mean(capacity_cluster_test_A);
    ergodic_cap_test_B(SNR_indx) = mean(capacity_cluster_test_B);
    ergodic_cap_test_C(SNR_indx) = mean(capacity_cluster_test_C);
    ergodic_cap_test_D(SNR_indx) = mean(capacity_cluster_test_D);
    
    
    ergodic_cap_test_E(SNR_indx) = mean(capacity_cluster_test_E);
    ergodic_cap_test_F(SNR_indx) = mean(capacity_cluster_test_F);
    ergodic_cap_test_G(SNR_indx) = mean(capacity_cluster_test_G);
    ergodic_cap_test_H(SNR_indx) = mean(capacity_cluster_test_H);
    
    ergodic_cap_test_I(SNR_indx) = mean(capacity_cluster_test_I);
    
    toc
end


figure(4)
hold on
p3 = plot(SNR_axis,ergodic_cap_test_A);
p3.LineWidth = 2;
p3 = plot(SNR_axis,ergodic_cap_test_B);
p3.LineWidth = 2;
p3 = plot(SNR_axis,ergodic_cap_test_C);
p3.LineWidth = 2;
p3 = plot(SNR_axis,ergodic_cap_test_D);
p3.LineWidth = 2;
hold off
legend('C = 20, L = 20', 'C = 20, L = 10', 'C = 20, L = 5', 'C = 20, L = 1')
ylabel('Ergodic Capcacity (bps/Hz)')
xlabel('SNR(dB)')
grid on


figure(5)
hold on
p3 = plot(SNR_axis,ergodic_cap_test_E);
p3.LineWidth = 2;
p3 = plot(SNR_axis,ergodic_cap_test_F);
p3.LineWidth = 2;
p3 = plot(SNR_axis,ergodic_cap_test_G);
p3.LineWidth = 2;
p3 = plot(SNR_axis,ergodic_cap_test_H);
p3.LineWidth = 2;
p3 = plot(SNR_axis,ergodic_cap_test_I);
p3.LineWidth = 2;
hold off
legend('C = 20, L = 20', 'C = 10, L = 20', 'C = 5, L = 20', 'C = 1, L = 20', ...
    'C = 1, L = 1')
ylabel('Ergodic Capcacity (bps/Hz)')
xlabel('SNR(dB)')
grid on



%-------------------------------------------



function y = rand_laplace(mu, sigma, m)    
    u = rand(m, 1)-0.5;
    b = sigma / sqrt(2);
    y = mu - b * sign(u).* log(1- 2* abs(u));
end













function capacity_classical = kronecker_based_model(SNR_dB, dim, rho)
    ant = 0:1:dim-1;
    ant_exponent = toeplitz(ant);
    phase_rho = abs(rho)*exp(1j*2*pi*rand(1,1)).*ones(dim,dim);
    Rx_corr = phase_rho.^ant_exponent;
    Rx_corr = Rx_corr / norm(Rx_corr,'fro'); % normalise correlation matrix

    % create random unit var, zero mean, normally distirbuted channel
    % normalise across 4 dimensional Gaussian to unit variance
    H_classical = sqrt(db2pow(SNR_dB))*sqrt(1/(2*dim))*(randn(dim) + 1j*randn(dim));
    
    term_classical = eye(dim) + H_classical*Rx_corr*ctranspose(H_classical);
    % add absolute term to make the logarithm behave
    capacity_classical = log2(abs(det(term_classical))); 
end








function capacity_cluster = cluster_based_model(rays_per_cluster, ...
    cluster_count, sigma_s, sigma_c, SNR_dB, dim, wavelength,  ...
    pathloss_cluster, cluster_dist_from_rx, cluster_shadow_var)

    antenna_dist= wavelength/2;

    
    
    cluster_shadowing = db2pow(cluster_shadow_var*randn(cluster_count,1));    
    cluster_dist = cluster_dist_from_rx * rand(cluster_count,1);
    
    beta_cluster = cluster_dist.^(-pathloss_cluster).*cluster_shadowing;
    
    
    
    
    
    %beta_cluster = 1;
    
    ant = 0:1:dim-1;

    summed_cluster_terms = zeros(dim,dim);
    
    
    cluster_coeff_phase = j*2*pi*rand(cluster_count,1); 
    
    gamma_c = sqrt(beta_cluster).*exp(cluster_coeff_phase)/cluster_count;
    
    gamma_c = gamma_c / (sum(gamma_c));
    
   

    for cluster_indx = 1:cluster_count
        

        
        
        central_cluster_angle = sigma_c*randn(1, 1);    


        summed_ray_terms = zeros(dim,dim);


        offset_angle = rand_laplace(0,sigma_s,rays_per_cluster);
        
        for ray_indx = 1:rays_per_cluster



            phi_AOD = central_cluster_angle + offset_angle(ray_indx);
            phi_AOA = central_cluster_angle + offset_angle(ray_indx);

            a_tx_AOD = exp(antenna_dist*j*2*pi*ant* cos(phi_AOD))';
            a_rx_AOA = exp(antenna_dist*j*2*pi*ant* cos(phi_AOA))';           

            h_iid_ray = (1\sqrt(2))*(randn(1,1) + j*randn(1,1));
            %h_iid_ray = 1;

            ray_term = h_iid_ray * a_rx_AOA * a_tx_AOD';
            
            
            summed_ray_terms = sqrt(gamma_c(cluster_indx) / rays_per_cluster) ...
             * ray_term + summed_ray_terms;



        end
        
        % normalise for cluster count
        summed_cluster_terms = summed_ray_terms + summed_cluster_terms;
        %summed_cluster_terms = summed_ray_terms + summed_cluster_terms;


    end
    summed_cluster_terms = summed_cluster_terms;
    
    H_cluster = sqrt(db2pow(SNR_dB))*(summed_cluster_terms);
       
    term = eye(dim) +  H_cluster*ctranspose(H_cluster);

    % add absolute term to make the logarithm behave    
    capacity_cluster = log2(abs(det(term))); 

end

























