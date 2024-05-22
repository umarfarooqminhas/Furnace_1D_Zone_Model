%% Initialization Function
function [Q_s_array,Q_g_array,Temp_s_array,Temp_mid_array,Temp_stock,Temp_stock_mid,eps_array,eta_fur,Temp_g_array] = Variables_Initialization(n_s,n_v,Bo_n)

    Q_s_array=zeros(n_s,length(Bo_n));
    Q_g_array=zeros(n_v,length(Bo_n));
    Temp_s_array=zeros(length(Bo_n),n_s);
    Temp_mid_array=zeros(length(Bo_n),n_s);
    Temp_stock=zeros(1,length(Bo_n));
    Temp_stock_mid=zeros(1,length(Bo_n));
    eps_array=zeros(1,length(Bo_n));
    eta_fur=zeros(1,length(Bo_n));
    Temp_g_array=zeros(length(Bo_n),n_v);

end