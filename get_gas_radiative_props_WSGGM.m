function [eps_g,abso,Beta,wi]=get_gas_radiative_props_WSGGM(vp,p,Temp_g)

    Lm=3.6*(vp(1,5)/p(1,8));                                                        % Mean beam length
    abso0=0;                                                                        %constant absorption coefficients
    abso1=0.267;
    abso2=4.65;
    abso3=72;
    a_g0=0.423+(0.0433*10^-3)*mean(Temp_g);                                         %Truelove temperature dependant weighting factors for WSGG model
    a_g1=0.285+(0.0513*10^-3)*mean(Temp_g);
    a_g2=0.227+(-0.0676*10^-3)*mean(Temp_g);
    a_g3=0.065+(-0.027*10^-3)*mean(Temp_g);
    eps_g=a_g0*(1-exp(-abso0*Lm))+a_g1*(1-exp(-abso1*Lm))+a_g2*(1-exp(-abso2*Lm))+a_g3*(1-exp(-abso3*Lm));  % WSGG model total gas emissivity
   
    eps_g=0.3;

    abso_g=(1/Lm)*log(1/(1-eps_g));                                                 % Effective gas absorption coefficient
   
    
   
    abso=abso_g;
    Beta = abso;                                                                    % Gas Volume Extinction coefficient
    wi=0/Beta;                                                                      % Non scattering Media - Scattering albedo =0

end