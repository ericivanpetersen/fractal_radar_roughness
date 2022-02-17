% Plot pretty plots!
close all;
clear all;
set(0,'defaultaxesfontsize',18);
set(0,'defaultlinelinewidth',3);

data                    = csvread('./Traces_Experiment/SHARAD_Raw3b.csv',1,1);
calib                   = 0;


% Anonymous function to calculate addition of decibels:
un_db = @(x) 10.^(x./10);
add_db = @(x,y) 10 .* log10( un_db(x) + un_db(y) );

% trim data?
%%std_cut                     = find(data(:,15)./data(:,14) < 0.5);
%data                        = data(std_cut,:);
point                        = data(:,6);
trim                         = find(point ~= 133 & point ~= 74);
data                         = data(trim,:);
% daycut                      = find(data(:,8) == 0);
% slopecut                    = find(data(:,10) < 2);
% data                        = data(slopecut,:);
%data                        = data(daycut,:);
%data                        = data(noenf,:);

%define data from pxy:
sub_pow                     = data(:,3);
surf_pow                    = data(:,4);
depth                       = data(:,5);
point                       = data(:,6);
site                        = data(:,7);
daytime                     = data(:,8);
rough_param                 = data(:,9);
slope                       = data(:,10);
sl_dat                      = data(:,11);
sl_std                      = data(:,12);
sl                          = data(:,13);
hurst_s                     = data(:,14);
hurst_r                     = data(:,15);
rough_param(find(rough_param == -99999)) = NaN;
sub_pow(find(sub_pow == -99999)) = NaN;
sub_pow_corr                = sub_pow + depth.*0.011;
sub_m_surf                  = sub_pow_corr - surf_pow;
NP                          = size(data,1);

r_eff                       = (5./(4.*pi^.2.*sl.^2.*cosd(0).^2)).^1./(2.*hurst_s);

%Site-specific indicies:
indA                       = find(site==1);
indB                       = find(site==2);
indD                       = find(site==4);
indC                       = find(site==3);
indE                       = find(site==5);
dind                       = [indA;indB;indC(find(isnan(sub_pow(indC) == 0)));indE]; %detection index
nind                       = [indD;indC(find(isnan(sub_pow(indC)) == 1))]; %non-detection at site D index
%calib_ind                  = find(point > 213 & point < 224);
calib_ind                   = find(point == 103);
if calib==1;
    calibration             = surf_pow(calib_ind) + 8.35
    surf_pow                = surf_pow - calibration;
    sub_pow                 = sub_pow - calibration;
    sub_pow_corr            = sub_pow_corr - calibration;
end

figure
plot(surf_pow(indA),r_eff(indA),'*',surf_pow(indB),r_eff(indB),'*',surf_pow(indC),r_eff(indC),'*',surf_pow(indD),r_eff(indD),'*',surf_pow(indE),r_eff(indE),'*');
legend('Site A','Site B','Site C','Site D','Site E','location','Southeast');
grid on
xlabel('Relative Surface Reflection Power, dB');
ylabel('r_{eff} (m)');
axis([-50 0 0 1000]);

figure
plot(sub_pow(indA),r_eff(indA),'*',sub_pow(indB),r_eff(indB),'*',sub_pow(indC),r_eff(indC),'*',sub_pow(indE),r_eff(indE),'*');
legend('Site A','Site B','Site C','Site D','Site E','location','Southeast');
grid on
xlabel('Relative Surface Reflection Power, dB');
ylabel('r_{eff} (m)');
axis([-50 0 0 1000]);

figure
plot(surf_pow(dind),r_eff(dind),'o','Markerfacecolor','b');
hold on
plot(surf_pow(nind),r_eff(nind),'rx','markersize',10);
legend('Detections','Nondetections','location','Southeast');
grid on
xlabel('Relative Surface Reflection Power (dB)');
ylabel('r_{eff} (m)');
axis([-45 -10 0 1000]);

figure
plot(surf_pow(indA),slope(indA),'*',surf_pow(indB),slope(indB),'*',surf_pow(indC),slope(indC),'*',surf_pow(indD),slope(indD),'*',surf_pow(indE),slope(indE),'*');
legend('Site A','Site B','Site C','Site D','Site E','location','Southeast');
grid on
xlabel('Relative Surface Reflection Power, dB');
ylabel('Slope, ^\circ');
axis([-50 0 0 6]);

figure
plot(surf_pow(dind),slope(dind),'o','Markerfacecolor','b');
hold on
plot(surf_pow(nind),slope(nind),'rx','markersize',10);
legend('Detections','Nondetections','location','Northeast');
grid on
xlabel('Relative Surface Reflection Power (dB)');
ylabel('Slope, ^\circ');
axis([-45 -10 0 6]);

%%
H                   = [0.1:0.01:0.6]';
s_l                 = [0.035:0.001:0.13];
rho                 = 0.1; %reflectivity corresponding to an eps = 5 surface.
theta               = 0;
theta_arr           = [0:0.1:8]';
TT                  = length(theta_arr);
HH                  = length(H);
SS                  = length(s_l);
sigma_0_map         = ones(HH,SS)*NaN;
sigma_0_map_slope   = ones(TT,SS)*NaN;
sigma_B_map         = sigma_0_map;
sigma_total_map     = sigma_0_map;

for ss=1:SS;
    for hh=1:HH;
        sigma_0_map(hh,ss)          = radar_backscatter_fractal_surface(rho,s_l(ss),H(hh),theta);
        sigma_B_map(hh,ss)          = calc_incoherent_backscatter_nadir_fractal(rho,s_l(ss),H(hh));
        sigma_total_map(hh,ss)      = add_db(sigma_0_map(hh,ss), sigma_B_map(hh,ss));
    end
    for tt=1:TT;
        sigma_0_map_slope(tt,ss)          = radar_backscatter_fractal_surface(rho,s_l(ss),0.4,theta_arr(tt));
    end
end

levels = [20,30,45,60,90,120,150,200]

map     = flipud(jet(100));

figure
pointsize                   = 50;
hold on
plot(atand(sl(indA)),hurst_s(indA),'*',atand(sl(indB)),hurst_s(indB),'*',atand(sl(indC)),hurst_s(indC),'*',atand(sl(indD)),hurst_s(indD),'*',atand(sl(indE)),hurst_s(indE),'*');
contour(atand(s_l),H,sigma_total_map,levels,'ShowText','off','LineColor','k')
xlabel('tan^{-1}s_\lambda (\circ)');
ylabel('H');
legend('Site A','Site B','Site C','Site D','Site E','location','Northeast');

figure
pointsize                   = 50;
%scatter(sl(dind),hurst_s(dind),pointsize,sub_pow_corr(dind),'filled')
hold on
plot(atand(sl(dind)),hurst_s(dind),'o', 'MarkerFaceColor', 'b');
plot(atand(sl(nind)),hurst_s(nind),'rx','MarkerSize',10);
contour(atand(s_l),H,sigma_total_map,levels,'ShowText','off','LineColor','k')
%colorbar
%colormap(map(30:100,:))
%caxis([min(sub_pow_corr(dind)) max(sub_pow_corr(dind))]);
xlabel('tan^{-1}s_\lambda (\circ)');
ylabel('H');
legend('Detections','Non-Detections','location','Southeast');

figure
scatter(atand(sl),hurst_s,pointsize,surf_pow,'filled');
colorbar
colormap(map)
hold on
contour(atand(s_l),H,sigma_total_map,levels,'ShowText','off','LineColor','k');
caxis([-40 max(surf_pow)]);
xlabel('tan^{-1}s_\lambda (\circ)');
ylabel('H');

%figure
%plot(sl(indA),sub_pow_corr(indA)-surf_pow(indA),'*',sl(indC),sub_pow_corr(indC) - surf_pow(indC),'*',sl(indE),sub_pow_corr(indE) - surf_pow(indE),'*');
%xlabel('RMS SLope');
%ylabel('Suburface Reflection Power');
%legend('Site A','Site C','Site E');

%%
eps2                            = 3; %Assumed dielectric constant of the debris layer              
theta_i_up                      = 2*slope - asind(1/sqrt(eps2)*sind(slope));
%theta_e_down                    = asind(sqrt(1/eps2*sind(slope)));
sl_scaled                       = sl./15.^hurst_s .* (15./sqrt(eps2)).^hurst_s;

sigma_0_s                       = ones(NP,1)*NaN;
sigma_0_slope0                  = sigma_0_s; 
sigma_b                         = sigma_0_s;
r_max                           = sigma_0_s;
rho                             = 0.1;
tau1                            = 1-rho;
tau2                            = tau1;
h                               = 300000;
D                               = 3000;
h_rms                           = sl * 15./sqrt(2);
for np=1:NP
    sigma_0_s(np)                       = radar_backscatter_fractal_surface(rho,sl(np),hurst_s(np),slope(np));
    sigma_0_slope0(np)                    = radar_backscatter_fractal_surface(rho,sl(np),hurst_s(np),0);
    [sigma_b(np), r_max(np)]                         = calc_incoherent_backscatter_nadir_fractal(rho,sl(np),hurst_s(np));
end

sigma_total = add_db(sigma_0_slope0, sigma_b)

coh_loss                   = 10*log10(exp(-2*(2*pi/15).^2*h_rms.^2.*cosd(slope).^2));
coh_loss_slope0            = 10*log10(exp(-2*(2*pi/15).^2*h_rms.^2.*cosd(0).^2));

r2_surf                         = corr(surf_pow,sigma_0_s).^2;
p_surf                          = polyfit(surf_pow,sigma_0_s,1);
[b_surf,bint_surf]              = regress(sigma_0_s-p_surf(2),surf_pow);
slope_err                       = diff(bint_surf);
x_range                         = [min(surf_pow)-10;max(surf_pow)+10];
y_model                         = polyval(p_surf,x_range);

r2_s0                         = corr(surf_pow,sigma_0_slope0).^2;
p_s0                          = polyfit(surf_pow,sigma_0_slope0,1);
y_mod0                         = polyval(p_s0,x_range);

r2_c                         = corr(surf_pow,coh_loss).^2;
p_c                          = polyfit(surf_pow,coh_loss,1);
[b_surfc,bint_surfc]              = regress(coh_loss-p_c(2),surf_pow);
y_modc                         = polyval(p_c,x_range);

r2_c0                         = corr(surf_pow,coh_loss_slope0).^2;
p_c0                          = polyfit(surf_pow,coh_loss_slope0,1);
y_modc0                         = polyval(p_c0,x_range);

%figure
%plot(surf_pow(indA),sigma_0_s(indA),'*',surf_pow(indB),sigma_0_s(indB),'*',surf_pow(indC),sigma_0_s(indC),'*',surf_pow(indD),sigma_0_s(indD),'*',surf_pow(indE),sigma_0_s(indE),'*');
%hold on
%plot(x_range,y_model,'k--');
%legend('Site A','Site B','Site C','Site D','Site E','location','Southeast');
%xlabel('Relative Surface Reflection Power (dB)');
%ylabel('\sigma_0 (dB)');
%grid on
%axis tight

figure
plot(surf_pow(indA),sigma_0_slope0(indA),'*',surf_pow(indB),sigma_0_slope0(indB),'*',surf_pow(indC),sigma_0_slope0(indC),'*',surf_pow(indD),sigma_0_slope0(indD),'*',surf_pow(indE),sigma_0_slope0(indE),'*');
hold on
%plot(x_range,y_mod0,'k--');
legend('Site A','Site B','Site C','Site D','Site E','location','Southeast');
xlabel('Relative Surface Reflection Power (dB)');
ylabel('\sigma_0 (dB)');
%title('Slope = 0');
grid on
axis([-50 0 20 160])

figure
scatter(surf_pow(dind),sigma_total(dind),'o','Markerfacecolor','b')
hold on
plot(surf_pow(nind),sigma_total(nind),'rx');
legend('Detections','Non-detections','location','Southeast');
xlabel('Relative Surface Reflection Power (dB)');
ylabel('\sigma_0 + \sigma_B (dB)');
%title('Slope = 0');
grid on
axis([-50 0 20 160])
%%
figure
scatter(surf_pow,sigma_total,pointsize,slope,'filled')
colorbar
colormap(flipud(map))
caxis([0, 4])
xlabel('Relative Surface Reflection Power (dB)');
ylabel('\sigma_0 + \sigma_B (dB)');
%title('Slope = 0');
grid on
axis([-50 0 20 160])

%%
figure
scatter(sigma_0_slope0(dind),sigma_b(dind),pointsize,'filled');
colorbar
colormap(flipud(map))
caxis([-40, -25])
hold on
plot(sigma_0_slope0(nind),sigma_b(nind),'rx','markersize',10);
legend('Detections','Non-detections','location','Southeast');
xlabel('\sigma_0 (dB)');
ylabel('\sigma_B (dB)');
%title('Slope = 0');
grid on

%%

figure
histogram(atand(sl(dind)),[2:0.25:7]);
hold on
histogram(atand(sl(nind)),[2:0.25:7]);
legend('Detections','Non-Detections')
xlabel('tan^{-1}s_\lambda (\circ)')
ylabel('Number of Observations')

figure
histogram(hurst_s(dind),[0.1:0.025:0.55]);
hold on
histogram(hurst_s(nind),[0.1:0.025:0.55]);
legend('Detections','Non-Detections')
xlabel('Hurst Exponent')
ylabel('Number of Observations')

figure
histogram(r_eff(dind),[50:50:1000]);
hold on
histogram(r_eff(nind),[50:50:1000]);
legend('Detections','Non-Detections')
xlabel('r_{eff} (m)')
ylabel('Number of Observations')

figure
histogram(sigma_0_slope0(dind),[30:5:160]);
hold on
histogram(sigma_0_slope0(nind),[30:5:160]);
legend('Detections','Non-Detections')
xlabel('\sigma_0 (dB)')
ylabel('Number of Observations')

figure
histogram(sigma_b(dind),[6:2:60]);
hold on
histogram(sigma_b(nind),[6:2:60]);
legend('Detections','Non-Detections')
xlabel('\sigma_B (dB)')
ylabel('Number of Observations')

figure
histogram(sigma_total(dind),[30:5:160]);
hold on
histogram(sigma_total(nind),[30:5:160]);
legend('Detections','Non-Detections')
xlabel('\sigma_0 + \sigma_B (dB)')
ylabel('Number of Observations')
grid on




