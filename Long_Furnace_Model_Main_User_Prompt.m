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

Furnace_Power = input("Enter Furnace overall power of all the burners in MW: "); 


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
Back=p(B_n,:);


%% Arrays Initialization
[Q_s_array,Q_g_array,Temp_s_array,Temp_mid_array,Temp_stock,Temp_stock_mid,eps_array,eta_fur,Temp_g_array] = Variables_Initialization(n_s,n_v,Bo_n);

%% Tranisent Loop

pos=Bo_n(1);
for w=1:length(Bo_n)
%% Zone Model Preliminary Calculations

    [TEAs,TEAs_SS,TEAs_SG,TEAs_GS,TEAs_GG,WWW,Q,P,C,R,X,D,bb,dd,ddd,eps_g,Beta] = Zone_Model_Preliminary_Calcs(n_v,n_s,n_t,vp,Temp_g,m,n,o,p,F_n,B_n,T_n,Bo_n,eps);

    %% Burners distribution inside the furnace
    Q_comb = burners_distrubution(Furnace_Power,N_burners,vp,n_v);                  % Burners distribution fuction in furnace volume zones
    
    %% Non-linear Energy Balance for Volume Zones
    [Temp_g] = Volume_Zone_Energy_Balance_Solver(TEAs_GG,TEAs_GS,Temp_s,vp,Beta,n_v,n_s,Total_flowrate,N_burners,Q_comb,Bo_n,p,B_n);
                                                                                    
    [Q_rad,Q_conv,Q_enth,Q_conv_1]=get_volume_zone_heat_transfer(TEAs_GG,TEAs_GS,Temp_g,Temp_s,vp,Bo_n,Total_flowrate,Q_comb,n_v,n_s,Beta,Back,N_burners); % Heat transfers from the volume zones calculations
                              
    %% Net Heat Transfer from Surface and Volume Zones
    [Q_s,Q_g]=Net_Heat_Transfer_Volume_Surface_Zones(eps,Temp_s,Temp_g',TEAs_SS,TEAs_SG,TEAs_GS,TEAs_GG,bb,dd);
    
    %% 1-D transient heat equation solution for all surface (walls/stock)
    [Temp_s,Temp_mid] = One_D_Transient_Heat_Conduction_base(Q_s,t_cycle,Bo_n,p,Temp_s,n_s);
        
    [eps_array,Temp_g_array,Temp_s_array,Temp_mid_array,eta_fur,Q_s_array,Q_g_array] = Store_Variables_Each_Instance(eps_g,Temp_g,Temp_s,Temp_mid,Q_rad,Q_comb,Q_s,Q_g,w,eps_array,Temp_g_array,Temp_s_array,Temp_mid_array,eta_fur,Q_s_array,Q_g_array);

    [Temp_stock,Temp_stock_mid,Temp_s,Temp_mid] = Update_Stock_Pos(Bo_n,Temp_s,Temp_mid,n_s,Temp_s_array,Temp_mid_array,w,pos,Temp_stock,Temp_stock_mid);
    
    disp(w);
    
    pos=pos+1;
end
