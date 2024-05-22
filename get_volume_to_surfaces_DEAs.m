  %% Volume to Surface DEA
  function VS = get_volume_to_surfaces_DEAs(m,n,o,vp,Beta,p,F_n,B_n,T_n,Bo_n)
 
    Front=p(F_n,:);
    Back=p(B_n,:);
    Bottom=p(Bo_n,:);
    Top=p(T_n,:);
    %Right=p(L_n,:);
    %Left=p(R_n,:);
    % DEAs from Volumes to all the other surfaces;
    VF=DEA_VS(Front,m*o,Beta,vp,m,n,o);
    VB=DEA_VS(Back,m*o,Beta,vp,m,n,o);
    VBo=DEA_VS(Bottom,n*m,Beta,vp,m,n,o);
    VT=DEA_VS(Top,n*m,Beta,vp,m,n,o);

    VS=[VF VB VBo VT];

  

function [VS]=DEA_VS(SS, s_nodes,Beta,vp,m,n,o)
node = 0;
VS=zeros(s_nodes,n*m*o);
    for p1=1:s_nodes % Number of nodes on the surface
        for vp1=1:n*m*o % Number of nodes in the enclosed volume
        node = node + 1;
        R= sqrt((vp(vp1,1)-SS(p1,1)).^2 + (vp(vp1,2)-SS(p1,2)).^2 + (vp(vp1,3)-SS(p1,3)).^2);
        cos1=(SS(p1,4)*sqrt((SS(p1,1)-vp(vp1,1)).^2)+SS(p1,5)*sqrt((SS(p1,2)-vp(vp1,2)).^2)+SS(p1,6)*sqrt((SS(p1,3)-vp(vp1,3)).^2))./R;
        VS(vp1,p1)=((exp(-(Beta).*R)).*((cos1.*Beta.*(SS(p1,7)).*(vp(vp1,4)))))./((R.^2).*(pi));
        end
    end
end

end