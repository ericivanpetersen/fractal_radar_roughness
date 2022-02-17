function [hurst,s_l,r2] = hurst_rmsslope(dx_array,slope_array,lambda,plot_num);
%Power-law fit to rms slope as a function of delta x to calculate the hurst
%exponent and rms slope at lambda.

%fit
hurst_fit                   = fit(dx_array,slope_array,'power1');
hurst                       = 1 + hurst_fit.b; %Hurst Exponent

%evaluate slope model:
slope_model             = hurst_fit.a .* dx_array.^hurst_fit.b;
s_l                     = slope_model(find(dx_array==lambda));

SSres = sum( (slope_array - slope_model).^2 );
SStot = sum( (slope_array - mean(slope_array)).^2 );

r2                      = 1 - SSres/SStot * (length(slope_array)-1)./(length(slope_array)-2);

%plot if desired:
if plot_num==1;
    figure
    loglog(dx_array,atand(slope_array),'*');
    hold on
    loglog(dx_array,atand(slope_model),'-');
    xlabel('\Deltax (m)');
    ylabel('tan^{-1}s (\circ)');
    legend('Data','Exponential Fit');
end    

end

