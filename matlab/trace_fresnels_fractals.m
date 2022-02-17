%Trace-by-trace experiment:

%Traces_Experiment All Surface!

%Get ready, kids!:
clear all;
close all;
addpath('../../../orig/');
addpath('./Traces_Experiment/RawDTM/');
set(0,'defaultaxesfontsize',20);
set(0,'defaultlinelinewidth',3);
outfile = fopen('./Traces_Experiment/SHARAD_Rough_out_1500m.csv','a');

%set bins, variables etc:
eps                         = [2:0.01:8];
powbin                      = [-50:0.1:-20];
bin                         = [0:0.1:15];
dx                          = [15:15:1500]';


%load in SHARAD Reflection Powers:
pxy                         = csvread('./Traces_Experiment/SHARAD_Raw2.csv',1,1);

%define data from pxy:
sub_pow                     = pxy(:,3);
surf_pow                    = pxy(:,4);
depth                       = pxy(:,5);
point                       = pxy(:,6);
site                        = pxy(:,7);
daytime                     = pxy(:,8);
rough_param                 = pxy(:,9);
slope                       = pxy(:,10);
rough_param(find(rough_param == -99999)) = NaN;
sub_pow(find(sub_pow == -99999)) = NaN;
sub_pow_corr                = sub_pow + depth.*0.0065;

%Site-specific indicies:
indA                       = find(site==1);
indB                       = find(site==2);
indD                       = find(site==4);
indC                       = find(site==3);
indE                       = find(site==5);
    
powA                        = surf_pow(find(site==1));
powB                        = surf_pow(find(site==2));
powD                        = surf_pow(find(site==4));
powC                        = surf_pow(find(site==3));
powE                        = surf_pow(find(site==5));

NP                          = size(pxy,1);
hurst_h                     = ones(NP,1)*NaN;
hurst_s                     = hurst_h;
sl                          = hurst_h;
sl_std                      = sl;
rms_h                       = sl;
rms_h_std                   = sl;
r2                          = sl;
sl_mod                      = sl;
%%

for np=1:NP;
    %np  = 15;
    file_str                = strcat('DTM',num2str(point(np),'%03.f'),'.tif');
    [sl(np),sl_std(np),sl_mod(np),hurst_s(np),r2(np)]   = DTM_rms_analy_2d(file_str,dx);
    fprintf('%s\r\n',file_str);
    fprintf('%4.5f,%4.5f,%4.5f,%4.5f,%4.5f\r\n',sl(np),sl_std(np),sl_mod(np),hurst_s(np),r2(np));
    fprintf(outfile,'%03.f,%4.5f,%4.5f,%4.5f,%4.5f,%4.5f\r\n',point(np),sl(np),sl_std(np),sl_mod(np),hurst_s(np),r2(np));
    close all
end

%%

detect_ind                  = find(isnan(sub_pow) == 0);
nondet_ind                  = indD;


figure
plot(sl_mod(detect_ind),hurst_s(detect_ind),'*',sl_mod(nondet_ind),hurst_s(nondet_ind),'*');
xlabel('RMS Slope');
ylabel('Hurst from RMS Slope');
legend('Detections','Non-Detections');

%%
sigma_0_h                       = ones(NP,1)*NaN;
sigma_0_h0                      = sigma_0_h;
sigma_0_s                       = sigma_0_h;
sigma_0_s_fore               = sigma_0_h;
sigma_0_h_fore               = sigma_0_h;
rho                             = 0.15;
for np=1:NP
    sigma_0_s(np)                       = radar_backscatter_fractal_surface(rho,sl_mod(np),hurst_s(np),slope(np));
    sigma_0_s_fore(np)                       = radar_scatter_fractal_surface(rho,sl_mod(np),hurst_s(np),slope(np),0);
end
sigma_0_s                       = sigma_0_s - median(sigma_0_s) + median(surf_pow);

coherent_loss                   = 10*log10(exp(-4*(2*pi/15*rms_h).^2));

figure
plot(surf_pow(indA),sigma_0_s(indA),'*',surf_pow(indB),sigma_0_s(indB),'*',surf_pow(indC),sigma_0_s(indC),'*',surf_pow(indD),sigma_0_s(indD),'*',surf_pow(indE),sigma_0_s(indE),'*');
legend('Site A','Site B','Site C','Site D','Site E');
xlabel('Surface Reflection Power');
ylabel('Radar Backscatter');

figure
plot(sub_pow_corr(indA),sigma_0_s_fore(indA),'*',sub_pow_corr(indB),sigma_0_s_fore(indB),'*',sub_pow_corr(indC),sigma_0_s_fore(indC),'*',sub_pow_corr(indD),sigma_0_s_fore(indD),'*',sub_pow_corr(indE),sigma_0_s_fore(indE),'*');

%%
site_surfpow_mean                = [mean(surf_pow(indA)); mean(surf_pow(indB)); ...
                                mean(surf_pow(indC)); mean(surf_pow(indD)); ...
                                mean(surf_pow(indE))];

site_surfpow_std                = [std(surf_pow(indA)); std(surf_pow(indB)); ...
                                std(surf_pow(indC)); std(surf_pow(indD)); ...
                                std(surf_pow(indE))];
                            
    
