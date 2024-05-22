clc;
clear;

%% Radiative heat transfer inside the furnace + Radiant section efficiency calcs. 
% Rectangular shaped furnace for steel slabs reheating

%% INPUTS Section
% Furnace Dimensions

L=input("Enter Furnace Length in meters: ");                                    % L = 36 Length in the x direction
W=input("Enter Furnace Width in meters: ");                                     % W = 10  Width in the y direction
H=input("Enter Furnace Height in meters: ");                                    % H = 2.5 Height in the direction


%% Fuel Inputs and Calculations

Furnace_Power = input("Enter Furnace overall power of all the burners in MW: "); % 20 MW
N_burners=input("Number of Burners: ");                                         % 20 Enter the number of burners present in furnace                                                   
t_cycle_minutes=input("Enter heating cycle time in minutes: ");                 % 90 mins Enter Furnace heating cycle time in minutes
t_cycle=t_cycle_minutes*60;
Fuel_comp = input("Please enter 1 for CH4 combustion and 2 for H2 combustion: ");
    while Fuel_comp ~= 1 && Fuel_comp ~= 2
        disp('Please enter a valid fuel input');
        Fuel_comp = input("Please enter 1 for CH4 combustion and 2 for H2 combustion: ");
    end
Oxidant_comp = input("Please enter 3 for combustion with Air and 4 for combustion with O2: ");
    while Oxidant_comp ~= 3 && Oxidant_comp ~= 4
        disp('Please enter a valid input');
        Oxidant_comp = input("Please enter 3 for combustion with Air and 4 for combustion with O2: ");
    end
excess = input("Please entire excess oxidant (Air or O2) percentage: ");
Total_flowrate = Fuel_calculations_per_burner(Furnace_Power,N_burners,Fuel_comp,Oxidant_comp,excess);  % Fuel + Oxidant flowrate calculation function per burner

%% Discretization of the Furnance
% Number of steps in the W, L and H direction:

n = 1;                                                                              % Increments in the W (y) direction
m = input("Enter number of Increments along the length of furnace (MUST be greater then nos. of burners): ");              % 90 Increments in the L (x) direction                                                                                  % Increments in the L (y) direction
o = 1;                                                                              % Increments in the H (z) direction
[n_t,n_s,n_v,F_n,B_n,Bo_n,T_n] = furnace_discretization_nodes(m,n,o);               % Furnace discretization into surface and volume zone function


%% %% FURNACE SURFACES + STOCK EMISIVITIES AND TEMPERATURE INITILIZATION

e_w=0.8;                                                                            % Walls emissivity considered as constant for simplicity 
eps=e_w*ones(1,n_s);                                                                % All surfaces emissivities
eps(1,Bo_n)=0.8;                                                                    % Stock emisivties considered as constant
Temp_s=1000*ones(n_s,1);                                                            % Surface Temperatures, top, front and back.(K)
Temp_s(Bo_n,1)=300;                                                                 % Stock temperatures (K)
Temp_g=300*ones(1,n_v);                                                             % Medium temperatures (K)


%% ZONE MODEL DEVELOPMENT

p = get_surface_nodes_properties(m,n,o,L,W,H);                                      % All surfaces nodes properties (x,y,z, normal vector, element area, entire surface area)
vp = get_volumes_nodes_properties(m,n,o,L,W,H);                                     % All Volume nodes properties



%% Arrays Initialization

Q_s_array=zeros(n_s,length(Bo_n));
Q_g_array=zeros(n_v,length(Bo_n));
Temp_s_array=zeros(length(Bo_n),n_s);
Temp_mid_array=zeros(length(Bo_n),n_s);
Temp_stock=zeros(1,length(Bo_n));
Temp_stock_mid=zeros(1,length(Bo_n));
eps_array=zeros(1,length(Bo_n));
eta_fur=zeros(1,length(Bo_n));
Temp_g_array=zeros(length(Bo_n),n_v);

%% Tranisent Heat Conduction Loop

pos=Bo_n(1);
for w=1:length(Bo_n)

    %% WEIGHTED SUM OF GRAY GASES MODEL
    [eps_g,abso,Beta,wi]=get_gas_radiative_props_WSGGM(vp,p,Temp_g);
    
    eps_array(w)=eps_g;   


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

    %% Burners distribution inside the furnace
    Q_comb = burners_distrubution(Furnace_Power,N_burners,vp,n_v);                  % Burners distribution fuction in furnace volume zones
        
    options = optimoptions(@fsolve,'Algorithm','trust-region','Display','iter','UseParallel',false,'OptimalityTolerance',1.0000);     % options for MATLAB fsolve non-linear equations solver
     
    initial_guess=zeros(1,n_v);                                                     % Initial guess for Temp_g, volumes zones temperature solution
    Back=p(B_n,:);
    [Temp_g] = fsolve(@(Temp_g)nonlinear_function_2(Temp_g,TEAs_GG,TEAs_GS,Temp_s,vp,Beta,n_v,n_s,Back,Total_flowrate,N_burners,Q_comb,Bo_n),initial_guess,options); % fsolve function call for non-linear energy equations solutions for Temp_g
                                              
    Temp_g_array(w,:)=Temp_g;                                                        % Volume Zones/ gas medium temperatures array at every instant of heating cycle

    [Q_rad,Q_conv,Q_enth,Q_conv_1]=get_volume_zone_heat_transfer(TEAs_GG,TEAs_GS,Temp_g,Temp_s,vp,Bo_n,Total_flowrate,Q_comb,n_v,n_s,Beta,Back,N_burners); % Heat transfers from the volume zones calculations

    eta_fur(w) = -sum(Q_rad)/sum(Q_comb);                                           % Overall radiant heat transfer to all the surfaces efficiency of radiant section
    
    Temp_g=Temp_g'; 
    
   
    Eb_s=(5.67*(10^-8))*(Temp_s.^4);                                                % Emissive powers of the surfaces [W/m2]
    Eb_g=(5.67*(10^-8))*(Temp_g.^4);                                                % Emissive powers of gas zones [W/m2]
    Q_s=(((eps.*bb)'.*Eb_s)-(TEAs_SS*Eb_s)-(TEAs_SG*Eb_g))/1000;                    % Heat trnasfer to all surfaces in kW
    Q_s_array(:,w)=Q_s;                                                             % Heat transfer to all surfaces array at each time instant
    Q_g=(dd'.*Eb_g-(TEAs_GS*Eb_s)-(TEAs_GG*Eb_g))/1000;                             % Heat trnasfer to all volumes zones in kW
    Q_g_array(:,w)=Q_g;                                                             % Heat transfer to/from all volumes array at each time instant
    
    %% 1-D transient heat equation solution for all surface (walls/stock)
    [Temp_s,Temp_mid] = One_D_Transient_Heat_Conduction_base(Q_s,t_cycle,Bo_n,p,Temp_s,n_s);
        
    Temp_s_array(w,:) = Temp_s;                                                       
     Temp_mid_array(w,:) = Temp_mid; 
    
    for i=1:n_s                                                                     % Update stock poistion loop by shifting the Temperature of previous surface node to next surface
         if i>Bo_n(1) && i<Bo_n(end)+1
            Temp_s(i)=Temp_s_array(w,i-1);
            Temp_mid(i)=Temp_mid_array(w,i-1);
        end
        if i==Bo_n(1)
             Temp_s(i)=300;                                                         % Stock temperature at initial position 
             Temp_mid(i)=300; 
        end
    end
    
    disp(w);
    Temp_stock(w)=Temp_s(pos);                                                      % Stock temperature array along the length of furnace for final plotting
    Temp_stock_mid(w)=Temp_mid(pos); 

    pos=pos+1;
end

