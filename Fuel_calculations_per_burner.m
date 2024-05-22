function Total_flowrate_per_burner = Fuel_calculations_per_burner(Furnace_Power,N_burners,Fuel_comp,Oxidant_comp,excess)

burner_pow=Furnace_Power/N_burners;
    if Fuel_comp == 1
        LHV_CH4 = 55;                                                               % LHV of CH4 in MJ/kg
        Fuel_flowrate_burner = burner_pow/LHV_CH4;
    end
    if Fuel_comp == 2
        LHV_H2 = 120;                                                               % LHV of H2 in MJ/kg
        Fuel_flowrate_burner = burner_pow/LHV_H2;
    end


    if Oxidant_comp == 3 && Fuel_comp == 1
            Oxidant_flowrate_burner = 1*2*(1+3.76)*Fuel_flowrate_burner/16*(1+excess/100)*28.84;
    end
    if Oxidant_comp == 3 && Fuel_comp == 2
            Oxidant_flowrate_burner = 1*0.5*(1+3.76)*Fuel_flowrate_burner/2*(1+excess/100)*28.84;
    end
    if Oxidant_comp == 4 && Fuel_comp == 1
            Oxidant_flowrate_burner = 1*2*Fuel_flowrate_burner/16*(1+excess/100)*32;
    end
    if Oxidant_comp == 4 && Fuel_comp == 2
            Oxidant_flowrate_burner = 1*0.5*Fuel_flowrate_burner/2*(1+excess/100)*32;
    end

Total_flowrate_per_burner = Fuel_flowrate_burner+Oxidant_flowrate_burner;                      %combined Fuel + Oxidant flowrate for simplified Heat transfer

end