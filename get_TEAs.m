function [TEAs,TEAs_SS,TEAs_SG,TEAs_GS,TEAs_GG,WWW,Q,P,C,R,X,D,bb,dd,ddd] = get_TEAs(p,Beta,abso,vp,eps,n_v,n_s,n_t,wi,DEA_ss,DEA_gg,DEA_gs,DEA_sg) 
bb=p(:,7)';                                                                   % Small areas of the surfaces of the discretised enclosure
    Area=ones(n_s,1)*bb;                                                            % Areas as a vector
    dd=4*Beta*vp(:,4)';                                                           % 4* extinction coef.*small volume
    ddd=4*abso*vp(:,4)';                                                          % 4* absorption coef.*small volume
    B=ones(n_v,1)*ddd;                                                              % 4* extinction coef.*small volume as a vector
    
    e=ones(n_s,1)*eps;                                                              % Emissivities as matrix
    ee=ones(n_v,1)*eps;
    
                                                  
    ww=wi*ones(1,n_v);                                                              % Scattering albedos
    wii=ones(n_v,1)*ww;
    www=ones(n_s,1)*ww;
    
   
    VF_ss=DEA_ss./Area;                                                              % Direct exchange areas to view factors conversion
    VF_gg=DEA_gg./B;
    
    T=zeros(n_s,n_s);
    % Developing the T matrix
    for i=1:n_s
        for j=1:n_s
            if j==i 
                T(i,j)=(1/e(i,j))-(VF_ss(i,j).*((1-e(i,j))./e(i,j)));
            else
                T(i,j)=(-1)*VF_ss(i,j).*((1-e(i,j))./e(i,j));
            end
        end
    end
    
    % S matrix
    S=DEA_ss.*e;
    
    % sg matrix
    %sg=DEA_sg;
    R=DEA_gs.*ee;
    Q=(DEA_gs./(ones(n_v,1)*bb)).*((1-ee)./ee);
    U=(DEA_sg./(ones(n_s,1)*ddd)).*www;
    WWW=zeros(n_v,n_v);
    for i=1:n_v
        for j=1:n_v
            if j==i 
            WWW(i,j)=(1/(1-wii(i,j)))-VF_gg(i,j).*wii(i,j);
            else
            WWW(i,j)=(-1)*VF_gg(i,j).*wii(i,j);
            end
        end
    end
    
    V=DEA_sg.*(1-www);
    X=DEA_gg.*(1-wii);
    C=S+U*(WWW\R);
    D=V+U*(WWW\X);
    P=T-U*(WWW\Q);
    
    % Total exchange areas
    YY=[T -U;-Q WWW];
    ZZ=[S V;R X];
    TEAs=(YY\ZZ);
    TEAs_SS=TEAs(1:n_s,1:n_s);
    TEAs_SG=TEAs(1:n_s,n_s+1:n_t);
    TEAs_GS=TEAs_SG';
    TEAs_GG=TEAs(n_s+1:n_t,n_s+1:n_t);
    %Test1=sum(TEAs(1,:));
    %Test2=sum(TEAs(:,n_t));
end