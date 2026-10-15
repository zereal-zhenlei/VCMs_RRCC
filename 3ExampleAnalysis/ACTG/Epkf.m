function [ y ] = Epkf( u  )
%EpKF Epanechnikov kernel function
y=3/4*(1-u.^2).*(abs(u)<=1);
%y=15/16*((1 - u.^2).^2).*(abs(u)<=1);
end

