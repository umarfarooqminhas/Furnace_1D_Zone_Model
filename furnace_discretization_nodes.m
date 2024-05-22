function [n_t,n_s,n_v,F_n,B_n,Bo_n,T_n] = furnace_discretization_nodes(m,n,o)



n_t=n*m*o+2*(n*o)+2*(m*o)+2*(n*m)-2;                                            % Total number of nodes (Surface + Volume), Left and right are entries and exit 
n_s=2*(n*o)+2*(m*o)+2*(n*m)-2;                                                  % Total number of surface nodes
n_v=n*m*o;                                                                      % Total number of volume nodes

%% surface node array
F_n=1:m*o;
B_n=m*o+1:2*(m*o);
Bo_n=2*(m*o)+1:2*(m*o)+n*m;
T_n=2*(m*o)+(n*m)+1:2*(m*o)+2*(n*m);

end
