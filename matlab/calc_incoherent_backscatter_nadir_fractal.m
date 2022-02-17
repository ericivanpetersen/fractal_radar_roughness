function [sigma_B, r_max] = calc_incoherent_backscatter_nadir_fractal(rho, s_l, H)
%Calculate the nadir incoherent backscattering coefficient for a fractal
%surface based on the models of Franceschetti et al (1999), Biccari et al.
%(2001) and Campbell and Shepard (2003)

sigma_B = rho./H .* ( (2*pi).^(H-1) ./ (s_l .* sqrt(2)) ).^(2/H) .* gamma(1./H);

r_max = sqrt( sigma_B ./ (4 * pi .* rho) );

sigma_B = 10 * log10(sigma_B);

end

