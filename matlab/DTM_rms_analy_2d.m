function [s_rms,s_rms_std,s_model,hurst_s,r2] = DTM_rms_analy_2d(DTM_file,dx)
%Version 2 cuts DTM profiles into 300 m segments

DTM                         = imread(DTM_file);
NaN_Const                   = DTM(1);
    
nind                        = find(DTM == NaN_Const);
DTM(nind)                   = NaN;
DTM(~any(~isnan(DTM), 2),:)     = [];
DTM(:,~any(~isnan(DTM), 1))     = [];

nan_map                     = isnan(DTM);

%% calculate rms slope for all profiles in DTM:

L                           = 15;
XX                          = length(dx);
%find all longitudinal/transvers profiles > 300 m
min_length                      = 300;
%min_length                      = max(dx);
prof_length_long                = size(nan_map,1) - sum(nan_map,1);
pind_lon                        = find(prof_length_long > min_length);
prof_length_trans               = size(nan_map,2) - sum(nan_map,2);
pind_trans                      = find(prof_length_trans > min_length);

PL                              = length(pind_lon);
PT                              = length(pind_trans);
sl_arr_y                        = ones(PL,XX)*NaN;
ad_arr_y                        = ones(PL,XX)*NaN;
sl_arr_x                        = ones(PT,XX)*NaN;
ad_arr_x                        = ones(PT,XX)*NaN;
h_rms_arr_y                     = ones(PL,XX)*NaN;
h_rms_arr_x                     = ones(PT,XX)*NaN;
H_ad_arr_x                      = ones(PT,1)*NaN;
H_ad_arr_y                      = ones(PT,1)*NaN;
for pl=1:PL;
    line                        = DTM(:,pind_lon(pl));
    prof                        = line(find(isnan(line)==0));
    [h_rms_arr_y(pl,:),ad_arr_y(pl,:),sl_arr_y(pl,:)] = profile_rough_analy([1:1:length(prof)]',prof,dx);
end

for pt=1:PT;
    line                        = DTM(pind_trans(pt),:);
    prof                        = line(find(isnan(line)==0));
    [h_rms_arr_x(pt,:),ad_arr_x(pt,:),sl_arr_x(pt,:)]   = profile_rough_analy([1:1:length(prof)]',prof',dx);
end

h_rms_arr                           = cat(1,h_rms_arr_x,h_rms_arr_y);
sl_arr                              = cat(1,sl_arr_x,sl_arr_y);
%isind                               = find(isnan(mean(h_rms_arr,2))==0);

mean_ad_arr                         = nanmean(cat(1,ad_arr_y,ad_arr_x));
std_ad_arr                          = nanstd(cat(1,ad_arr_y,ad_arr_x));
mean_sl_arr                         = nanmean(cat(1,sl_arr_y,sl_arr_x));
std_sl_arr                          = nanstd(cat(1,sl_arr_y,sl_arr_x));
h_rms_mean                          = nanmean(cat(1,h_rms_arr_x,h_rms_arr_y));
h_rms_std_arr                       = nanstd(cat(1,h_rms_arr_x,h_rms_arr_y));

xind                                = find(dx == L);
h_rms                               = h_rms_mean(xind);
h_rms_std                           = h_rms_std_arr(xind);
s_rms                               = mean_sl_arr(xind);
s_rms_std                           = std_sl_arr(xind);

[hurst_s,s_model,r_s]                     = hurst_rmsslope(dx,mean_sl_arr',15,1);

r2                                  = r_s;

end

