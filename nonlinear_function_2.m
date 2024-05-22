%% Function for solving energy balance of all volume nodes using MATLAB fsolve

function ObjF=nonlinear_function_2(Temp_g,TEAs_GG,TEAs_GS,Temp_s,vp,Beta,n_v,n_s,Back,Total_flowrate,N_burners,Q_comb,Bo_n)

%% Radiant heat transfer from each zone
Q_rad_g=zeros(1,n_v);
Q_rad_s=zeros(1,n_v);
Q_rad=zeros(1,n_v);
for i=1:n_v
    for j=1:n_v
    Q_rad_g(j)= 5.67*(10^-8)*TEAs_GG(i,j)*Temp_g(1,j)^4;                               % Radiant Heat transfer between gas zones
    end

     for j=1:n_s
    Q_rad_s(j)= 5.67*(10^-8)*TEAs_GS(i,j)*Temp_s(j,1)^4;                               % Radiant Heat transfer between gas and suface zones
     end
    Q_rad(i)=sum(Q_rad_g)+sum(Q_rad_s)-4*vp(i,4)*Beta*5.67*(10^-8)*Temp_g(1,i)^4;
end

%% Convective heat transfer from each zone
Q_conv=zeros(1,n_v);
h=5.6;
A=Back(1,7); %%CHECK PLEASE
for i=1:n_v                                                                             % Convective Heat transfer considered only with top and bottom surface of furnace for simplicity
        Q_conv_1=h*A*(Temp_g(1,i)-Temp_s(i+Bo_n(1)-1,1));
        Q_conv_2=h*A*(Temp_g(1,i)-Temp_s(i+Bo_n(end),1));
        Q_conv(i)=Q_conv_1+Q_conv_2;
end

%% Combustion Heat transfer at the zones where buners are installed

%Q_comb;                                       % Dividing total fuance power by each burner and considering uniform heat source at entire volume zone where burner is installed

%% Net Enthalpy flow from each zone
m_gas=Total_flowrate; 
cp_gas=1020;
   %% Initializations
burner_pos=zeros(1,N_burners);
Q_enth_out=zeros(1,n_v);
Q_enth_in=zeros(1,n_v);
Q_enth=zeros(1,n_v);


j=1;
for i=1:n_v
    if Q_comb(i)~=0
        burner_pos(j)=i;
        j=j+1;
    end
end
   Q_enth(1)=-m_gas*cp_gas*(Temp_g(1));                                                 % Net enthalpy flow at first zone, incoming fuel and air heating not considered

for i=1:n_v
    if i>1
        Q_enth_in(i)= m_gas*cp_gas*(Temp_g(i-1)); 
    
    for burner_check=1:length(burner_pos)
        if i>1
            if i==burner_pos(burner_check)
                m_gas=m_gas+Total_flowrate; 
            end
        end
    end
    
        Q_enth_out(i)= m_gas*cp_gas*(-Temp_g(i));                                        % Net enthalpy flow at last zone.
    
    Q_enth(i)=Q_enth_in(i)+Q_enth_out(i);
    end
end

ObjF = Q_comb+Q_enth+Q_rad-Q_conv;
end