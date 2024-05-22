function [Q_s,Q_g]=Net_Heat_Transfer_Volume_Surface_Zones(eps,Temp_s,Temp_g,TEAs_SS,TEAs_SG,TEAs_GS,TEAs_GG,bb,dd)
    Eb_s=(5.67*(10^-8))*(Temp_s.^4);                                                % Emissive powers of the surfaces [W/m2]
    Eb_g=(5.67*(10^-8))*(Temp_g.^4);                                                % Emissive powers of gas zones [W/m2]
    Q_s=(((eps.*bb)'.*Eb_s)-(TEAs_SS*Eb_s)-(TEAs_SG*Eb_g))/1000;                    % Heat trnasfer to all surfaces in kW
    Q_g=(dd'.*Eb_g-(TEAs_GS*Eb_s)-(TEAs_GG*Eb_g))/1000;                             % Heat trnasfer to all volumes zones in kW

end
