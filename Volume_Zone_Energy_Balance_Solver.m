function [Temp_g] = Volume_Zone_Energy_Balance_Solver(TEAs_GG,TEAs_GS,Temp_s,vp,Beta,n_v,n_s,Total_flowrate,N_burners,Q_comb,Bo_n,p,B_n)

options = optimoptions(@fsolve,'Algorithm','trust-region','Display','iter','UseParallel',false,'OptimalityTolerance',1.0000);     % options for MATLAB fsolve non-linear equations solver
    
    initial_guess=zeros(1,n_v);                                                                                                     % Initial guess for Temp_g, volumes zones temperature solution
    Back=p(B_n,:);
    Temp_g = fsolve(@(Temp_g)nonlinear_function_2(Temp_g,TEAs_GG,TEAs_GS,Temp_s,vp,Beta,n_v,n_s,Back,Total_flowrate,N_burners,Q_comb,Bo_n),initial_guess,options); % fsolve function call for non-linear energy equations solutions for Temp_g
                                              
    
end
