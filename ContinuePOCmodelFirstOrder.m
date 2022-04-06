clearvars
close all

load('POC model First Order Outputs.mat')

% k_sorp = 0.7/365;                 %Lerner
% k_desorp = 2/365;                 %Lerner
starttimes = [9.5,50];
endtimes = [40,137];
% jmp_k_sorp = k_sorp/5;
% jmp_k_desorp = k_desorp/5;
numiter=500*10^3;

% %Calculating values for log-normal priors for k_sorp (assuming Coefficient of Variation = 0.5)
% CV=0.5;
% expected_mu_k_sorp = k_sorp;
% var_mu_k_sorp = k_sorp^2*CV^2;
% mu_k_sorp=log(expected_mu_k_sorp^2/sqrt(var_mu_k_sorp+expected_mu_k_sorp^2));
% sig_k_sorp=sqrt(log(var_mu_k_sorp/expected_mu_k_sorp^2+1));
% 
% %Calculating values for log-normal priors for k_desorp (assuming Coefficient of Variation = 0.5)
% CV=0.5;
% expected_mu_k_desorp = k_desorp;
% var_mu_k_desorp = k_desorp^2*CV^2;
% mu_k_desorp=log(expected_mu_k_desorp^2/sqrt(var_mu_k_desorp+expected_mu_k_desorp^2));
% sig_k_desorp=sqrt(log(var_mu_k_desorp/expected_mu_k_desorp^2+1));
% 
% % x=0:0.00001:k_sorp*5
% % y=lognpdf(x,mu_k_sorp,sig_k_sorp)
% % plot(x,y)
% % hold on
% % plot([k_sorp,k_sorp],[0 100],'r')
% 
% %Initial (single) model run to find log likelihood and prior
% [norm_misfit0] = POCmodelFirstOrder(k_sorp,k_desorp,starttimes,endtimes,0);
% 
% loglike0 = -1/2 * sum(norm_misfit0.^2);
% prior0 = lognpdf(k_sorp,mu_k_sorp,sig_k_sorp)*lognpdf(k_desorp,mu_k_desorp,sig_k_desorp);
% prob0 = exp(loglike0)*prior0;
% 
% 
% %Markov Chain Monte Carlo
% accepted=0;
% rejected=0;
% k_sorp_track = k_sorp;
% k_desorp_track = k_desorp;
% loglike_track = loglike0;
% prior_track = prior0;
for i=1:numiter
    i
    tic
    proposal_k_sorp = k_sorp + randn*jmp_k_sorp;
    proposal_k_desorp = k_desorp + randn*jmp_k_desorp;
    
    [proposal_norm_misfit] = POCmodelFirstOrder(proposal_k_sorp,proposal_k_desorp,starttimes,endtimes,0);
    
    proposal_loglike = -1/2 * sum(proposal_norm_misfit.^2);
    proposal_prior = lognpdf(proposal_k_sorp,mu_k_sorp,sig_k_sorp)*lognpdf(proposal_k_desorp,mu_k_desorp,sig_k_desorp);
    
    prob = exp(proposal_loglike - loglike0)*proposal_prior/prior0;
    
    if rand<prob
        
        accepted = accepted+1;
        k_sorp = proposal_k_sorp;
        k_desorp = proposal_k_desorp;
        loglike0 = proposal_loglike;
        norm_misfit0 = proposal_norm_misfit;
        prior0 = proposal_prior;
    else
        rejected = rejected+1;
    end
    k_sorp_track(end+1) = k_sorp;
    k_desorp_track(end+1) = k_desorp;
    loglike_track(end+1) = loglike0;
    prior_track(end+1) = prior0;
    
    if mod(i,1000)==0
        origin = 'RunPOCmodelFirstOrder.m'
        save('POC model First Order Outputs.mat','origin','k_sorp_track','k_desorp_track','loglike_track','prior_track',...
            'accepted','rejected','k_sorp','k_desorp','loglike0','norm_misfit0','prior0',...
            'jmp_k_sorp','jmp_k_desorp','mu_k_sorp','mu_k_desorp','sig_k_sorp','sig_k_desorp')
    end
    time1=toc
end

acceptedratio = accepted/(accepted+rejected);

figure('Position',[50 50 1000 600])
subplot(2,2,1)
plot(loglike_track)
title('Log Likelihood')
subplot(2,2,2)
plot(k_sorp_track)
title('k_s_o_r_p')
ylabel('k_s_o_r_p (d^-^1)')
subplot(2,2,3)
plot(k_desorp_track)
title('k_d_e_s_o_r_p')
ylabel('k_d_e_s_o_r_p (d^-^1)')