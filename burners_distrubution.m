function Q_comb = burners_distrubution(Furnace_Power,N_burners,vp,n_v)
HeatSource=Furnace_Power/N_burners/vp(1,4)'*1E6;
     Q_comb=zeros(1,n_v);
        i=1;

    while i<n_v
        Q_comb(1,i)=HeatSource*vp(1,4)';
        i=i+round(n_v/N_burners);
    end

end