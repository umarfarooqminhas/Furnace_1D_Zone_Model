function [Temp_s,Temp_mid] = One_D_Transient_Heat_Conduction_base(Q_s,t_cycle,Bo_n,p,Temp_s,n_s)
Temp_mid=Temp_s;                                                                    %Initialization
    for i=1:n_s                                                                     % All surfaces loop
        if i<Bo_n(1)                                                                % Front back surface loop                                                                            
           fa=0;                                                                    % Adiabatic condition at furnace outer wall for simplicity
           fb=-Q_s(i)*1000/p(i,7);                                                  % Heat flux on furnae inner walls
           Tin=Temp_s(i);                                                           % Initial temepratures
           t=t_cycle/length(Bo_n);                                                  % time delay considered at each position, considered same heat flux for that duration
           T4_FB = One_D_Trasient_Heat_Conduction_CN_Refractory(fa,fb,Tin,t);       % Transient Heat equation solver function call to get end time temperatures
           Temp_s(i)=T4_FB(end,end);                                                % Update temperatures for next time step
           %{
           if Temp_s(i)>1800
               Temp_s(i)=1800;
           end
           %}
           Temp_mid(i)=T4_FB(round(end/2),round(end/2)); 
        end

        if i>Bo_n(1)-1 && i<Bo_n(end)+1                                             % stock temperature calculation loop
                 fa=-Q_s(i)*1000/p(i,7);                                            % Upper sider heat flux
               fb=-Q_s(i)*1000/p(i,7);                                              % Lower side heat flux (symmetry considered for simplicity
               Tin=Temp_s(i);
               t=t_cycle/length(Bo_n);
            T4_S = One_D_Trasient_Heat_Conduction_CN(fa,fb,Tin,t);
            Temp_s(i)=T4_S(end,end); 
            Temp_mid(i)=T4_S(round(end/2),round(end/2)); 
        end

          if i>Bo_n(end) && i<n_s(end)+1                                            % Top surface loop                                 
                 fa=0;
               fb=-Q_s(i)*1000/p(i,7);
               Tin=Temp_s(i);
               t=t_cycle/length(Bo_n);
            T4_T = One_D_Trasient_Heat_Conduction_CN_Refractory(fa,fb,Tin,t);
            Temp_s(i)=T4_T(end,end); 
           %{
           if Temp_s(i)>1800
               Temp_s(i)=1800;
           end
           %}
            Temp_mid(i)=T4_T(round(end/2),round(end/2)); 
          end
    end
end
