function [eps_array,Temp_g_array,Temp_s_array,Temp_mid_array,eta_fur,Q_s_array,Q_g_array] = Store_Variables_Each_Instance(eps_g,Temp_g,Temp_s,Temp_mid,Q_rad,Q_comb,Q_s,Q_g,w,eps_array,Temp_g_array,Temp_s_array,Temp_mid_array,eta_fur,Q_s_array,Q_g_array)
    eps_array(w)=eps_g;     
    Temp_g_array(w,:)=Temp_g;    
    Temp_s_array(w,:) = Temp_s;                                                       
    Temp_mid_array(w,:) = Temp_mid;
    eta_fur(w) = -sum(Q_rad)/sum(Q_comb); 
    Q_s_array(:,w)=Q_s;                                                             % Heat transfer to all surfaces array at each time instant
    Q_g_array(:,w)=Q_g;                                                             % Heat transfer to/from all volumes array at each time instant
    
end
