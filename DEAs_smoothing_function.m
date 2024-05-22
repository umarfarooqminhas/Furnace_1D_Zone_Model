function XXX=DEAs_smoothing_function(SS,VS,VV,Beta,p,vp,n_t)
                                                                                    % LEAST SQUARES SMOOTHING USING LAGRANGE
    XX=[SS VS';VS VV];                                                              % DEAs matrix
    W=XX.^2;                                                                        % Weights Matrix hat allow for DEAs to be adjusted proportionally
    AA=zeros(1,n_t);
    WW=zeros(1,n_t);
                                                                                    % Totals along the rows and columns
    for i=1:n_t
        for j=1:n_t
            AA(1,j)=sum(XX(:,j));                                                    % Sum of the rows of every column of the DEAs matrix
            WW(1,j)=sum(W(:,j));                                                     % Sum of the rows of every column of the weights matrix
        end
    end
                                                                                    % The conservation constraints that must be satisfied
    b=[p(:,7)' (4*Beta*vp(:,4))'];
    C=b-AA;
                                                                                    % Weights used to calculate the Lagrange multipliers
    R=W+(ones(n_t,1)*WW).*eye(n_t,n_t);
    
    L=R\C';                                                                         % Lagrange multipliers
 
    L1=ones(n_t,1)*L';                                                              % Lagrangian transpose for every DEA (rows)
  
    L2=L*(ones(1,n_t));                                                             % Lagrangian for every DEA (columns)
    
    lam=L1+L2;                                                                      % Lagrange multipliers for every DEA
    
    XXX=XX+W.*lam;                                                                  % Least squares smoothing using Lagrange multipliers (Smoothed DEAs)

end
    
  
    
    