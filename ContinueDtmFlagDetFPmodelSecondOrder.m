clearvars
close all

load('Dtm-Flag-Det-FP model Second Order Outputs.mat')

% k_sorp_dtm = 0.013;                 %m3 / mmol C / d, Stukel et al. (2019)
% k_sorp_flag = 0.013;                 %m3 / mmol C / d, Stukel et al. (2019)
% k_sorp_det = 0.013;                 %m3 / mmol C / d, Stukel et al. (2019)
% k_desorp = 2/365;                 %Lerner
% Eg_Th = 0.5;
starttimes = [9.5,50];
endtimes = [40,137];
% jmp_k_sorp_dtm = k_sorp_dtm/30;
% jmp_k_sorp_flag = k_sorp_flag/30;
% jmp_k_sorp_det = k_sorp_det/30;
% jmp_k_desorp = k_desorp/10;
% jmp_Eg_Th = 0.005;
numiter=1000*10^3;

% %Calculating values for log-normal priors for k_sorp_dtm (assuming Coefficient of Variation = 0.5)
% CV=0.5;
% expected_mu_k_sorp_dtm = k_sorp_dtm;
% var_mu_k_sorp_dtm = k_sorp_dtm^2*CV^2;
% mu_k_sorp_dtm=log(expected_mu_k_sorp_dtm^2/sqrt(var_mu_k_sorp_dtm+expected_mu_k_sorp_dtm^2));
% sig_k_sorp_dtm=sqrt(log(var_mu_k_sorp_dtm/expected_mu_k_sorp_dtm^2+1));
% 
% %Calculating values for log-normal priors for k_sorp_flag (assuming Coefficient of Variation = 0.5)
% CV=0.5;
% expected_mu_k_sorp_flag = k_sorp_flag;
% var_mu_k_sorp_flag = k_sorp_flag^2*CV^2;
% mu_k_sorp_flag=log(expected_mu_k_sorp_flag^2/sqrt(var_mu_k_sorp_flag+expected_mu_k_sorp_flag^2));
% sig_k_sorp_flag=sqrt(log(var_mu_k_sorp_flag/expected_mu_k_sorp_flag^2+1));
% 
% %Calculating values for log-normal priors for k_sorp_det (assuming Coefficient of Variation = 0.5)
% CV=0.5;
% expected_mu_k_sorp_det = k_sorp_det;
% var_mu_k_sorp_det = k_sorp_det^2*CV^2;
% mu_k_sorp_det=log(expected_mu_k_sorp_det^2/sqrt(var_mu_k_sorp_det+expected_mu_k_sorp_det^2));
% sig_k_sorp_det=sqrt(log(var_mu_k_sorp_det/expected_mu_k_sorp_det^2+1));
% 
% %Calculating values for log-normal priors for k_desorp (assuming Coefficient of Variation = 0.5)
% CV=0.5;
% expected_mu_k_desorp = k_desorp;
% var_mu_k_desorp = k_desorp^2*CV^2;
% mu_k_desorp=log(expected_mu_k_desorp^2/sqrt(var_mu_k_desorp+expected_mu_k_desorp^2));
% sig_k_desorp=sqrt(log(var_mu_k_desorp/expected_mu_k_desorp^2+1));
% 
% %Calculating values for beta distriubtion priors for Eg_Th (assuming Coefficient of Variation = 0.5)
% CV = 0.5;
% expected_Eg_Th = Eg_Th;
% var_Eg_Th = Eg_Th^2*CV^2;
% alpha_Eg_Th = (expected_Eg_Th*(1-expected_Eg_Th)/var_Eg_Th - 1) * expected_Eg_Th;
% beta_Eg_Th =  (expected_Eg_Th*(1-expected_Eg_Th)/var_Eg_Th - 1) * (1-expected_Eg_Th);
% 
% % x=0:0.00001:k_sorp*5
% % y=lognpdf(x,mu_k_sorp,sig_k_sorp)
% % plot(x,y)
% % hold on
% % plot([k_sorp,k_sorp],[0 100],'r')
% 
% % x = 0:0.001:1
% % y = betapdf(x,alpha_Eg_Th,beta_Eg_Th)
% % plot(x,y)
% % hold on
% % plot([Eg_Th,Eg_Th],[1 1],'r')
% 
% 
% %Initial (single) model run to find log likelihood and prior
% [norm_misfit0] = DtmFlagDetFPmodelSecondOrder(k_sorp_dtm,k_sorp_flag,k_sorp_det,k_desorp,Eg_Th,starttimes,endtimes,0);
% 
% loglike0 = -1/2 * sum(norm_misfit0.^2);
% prior0 = lognpdf(k_sorp_dtm,mu_k_sorp_dtm,sig_k_sorp_dtm)*lognpdf(k_sorp_flag,mu_k_sorp_flag,sig_k_sorp_flag)*lognpdf(k_sorp_det,mu_k_sorp_det,sig_k_sorp_det)*lognpdf(k_desorp,mu_k_desorp,sig_k_desorp)* betapdf(Eg_Th,alpha_Eg_Th,beta_Eg_Th);
% prob0 = exp(loglike0)*prior0;
% 
% 
% %Markov Chain Monte Carlo
% accepted=0;
% rejected=0;
% k_sorp_dtm_track = k_sorp_dtm;
% k_sorp_flag_track = k_sorp_flag;
% k_sorp_det_track = k_sorp_det;
% k_desorp_track = k_desorp;
% Eg_Th_track = Eg_Th;
% loglike_track = loglike0;
% prior_track = prior0;
for i=1:numiter
    i
    tic
    proposal_k_sorp_dtm = k_sorp_dtm + randn*jmp_k_sorp_dtm;
    proposal_k_sorp_flag = k_sorp_flag + randn*jmp_k_sorp_flag;
    proposal_k_sorp_det = k_sorp_det + randn*jmp_k_sorp_det;
    proposal_k_desorp = k_desorp + randn*jmp_k_desorp;
    proposal_Eg_Th = Eg_Th + randn*jmp_Eg_Th;
    
    [proposal_norm_misfit] = DtmFlagDetFPmodelSecondOrder(proposal_k_sorp_dtm,proposal_k_sorp_flag,proposal_k_sorp_det,proposal_k_desorp,proposal_Eg_Th,starttimes,endtimes,0);
    
    proposal_loglike = -1/2 * sum(proposal_norm_misfit.^2);
    proposal_prior = lognpdf(proposal_k_sorp_dtm,mu_k_sorp_dtm,sig_k_sorp_dtm)*lognpdf(proposal_k_sorp_flag,mu_k_sorp_flag,sig_k_sorp_flag)*lognpdf(proposal_k_sorp_det,mu_k_sorp_det,sig_k_sorp_det)*lognpdf(proposal_k_desorp,mu_k_desorp,sig_k_desorp)*betapdf(proposal_Eg_Th,alpha_Eg_Th,beta_Eg_Th);
    
    prob = exp(proposal_loglike - loglike0)*proposal_prior/prior0;
    
    if rand<prob
        
        accepted = accepted+1;
        k_sorp_dtm = proposal_k_sorp_dtm;
        k_sorp_flag = proposal_k_sorp_flag;
        k_sorp_det= proposal_k_sorp_det;
        k_desorp = proposal_k_desorp;
        Eg_Th = proposal_Eg_Th;
        loglike0 = proposal_loglike;
        norm_misfit0 = proposal_norm_misfit;
        prior0 = proposal_prior;
    else
        rejected = rejected+1;
    end
    k_sorp_dtm_track(end+1) = k_sorp_dtm;
    k_sorp_flag_track(end+1) = k_sorp_flag;
    k_sorp_det_track(end+1) = k_sorp_det;
    k_desorp_track(end+1) = k_desorp;
    Eg_Th_track(end+1) = Eg_Th;
    loglike_track(end+1) = loglike0;
    prior_track(end+1) = prior0;
    
    if mod(i,100)==0
        origin = 'RunDtmFlagDetFPmodelSecondOrder.m'
        save('Dtm-Flag-Det-FP model Second Order Outputs.mat','origin','k_sorp_dtm_track','k_sorp_flag_track','k_sorp_det_track','k_desorp_track','Eg_Th_track','loglike_track','prior_track',...
            'accepted','rejected','k_sorp_dtm','k_sorp_flag','k_sorp_det','k_desorp','Eg_Th','loglike0','norm_misfit0','prior0',...
            'jmp_k_sorp_dtm','jmp_k_sorp_flag','jmp_k_sorp_det','jmp_k_desorp','jmp_Eg_Th','mu_k_sorp_dtm','mu_k_sorp_flag','mu_k_sorp_det','mu_k_desorp','sig_k_sorp_dtm','sig_k_sorp_flag','sig_k_sorp_det','sig_k_desorp','alpha_Eg_Th','beta_Eg_Th')
    end
    time1=toc
end

acceptedratio = accepted/(accepted+rejected);

figure('Position',[50 50 1000 600])
subplot(2,3,1)
plot(loglike_track)
title('Log Likelihood')
subplot(2,3,2)
plot(k_sorp_dtm_track)
title('k_s_o_r_p_,_d_t_m')
ylabel('k_s_o_r_p (d^-^1)')
subplot(2,3,3)
plot(k_sorp_flag_track)
title('k_s_o_r_p_,_f_l_a_g')
ylabel('k_s_o_r_p (d^-^1)')
subplot(2,3,4)
plot(k_sorp_det_track)
title('k_s_o_r_p_,_d_e_t')
ylabel('k_s_o_r_p (d^-^1)')
subplot(2,3,5)
plot(k_desorp_track)
title('k_d_e_s_o_r_p')
ylabel('k_d_e_s_o_r_p (d^-^1)')
subplot(2,3,6)
plot(k_desorp_track)
title('Eg_T_h_,_t_r_a_c_k')
ylabel('Eg_T_h_,_t_r_a_c_k')