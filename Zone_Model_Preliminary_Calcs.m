function [TEAs,TEAs_SS,TEAs_SG,TEAs_GS,TEAs_GG,WWW,Q,P,C,R,X,D,bb,dd,ddd,eps_g,Beta] = Zone_Model_Preliminary_Calcs(n_v,n_s,n_t,vp,Temp_g,m,n,o,p,F_n,B_n,T_n,Bo_n,eps)
%% WEIGHTED SUM OF GRAY GASES MODEL
    [eps_g,abso,Beta,wi]=get_gas_radiative_props_WSGGM(vp,p,Temp_g);
    
     


    %% DIRECT EXCHANGE AEREAS CALCULATION
    SS = get_surfaces_DEAs(p,m,n,o,Beta,F_n,B_n,T_n,Bo_n);
    VV = get_volume_DEAs(vp,m,n,o,Beta);
    VS = get_volume_to_surfaces_DEAs(m,n,o,vp,Beta,p,F_n,B_n,T_n,Bo_n);

 
    
    %% Smoothing of DEAs
    XXX=DEAs_smoothing_function(SS,VS,VV,Beta,p,vp,n_t);

  
    DEA_ss=XXX(1:n_s,1:n_s);                                                        % Surface to surface DEAs
    DEA_sg=XXX(1:n_s,n_s+1:n_t);                                                    % Volume to surface DEAs
    DEA_gs=XXX(n_s+1:n_t,1:n_s);                                                    % Volume to surface DEAs
    DEA_gg=XXX(n_s+1:n_t,n_s+1:n_t);                                                % Volume to volume DEAs
    
    
    
    %% Total Exchange Area calculations
    [TEAs,TEAs_SS,TEAs_SG,TEAs_GS,TEAs_GG,WWW,Q,P,C,R,X,D,bb,dd,ddd] = get_TEAs(p,Beta,abso,vp,eps,n_v,n_s,n_t,wi,DEA_ss,DEA_gg,DEA_gs,DEA_sg);
end
