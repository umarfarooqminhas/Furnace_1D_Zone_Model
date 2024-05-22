function SS = get_surfaces_DEAs(p,m,n,o,Beta,F_n,B_n,T_n,Bo_n)
    Front=p(F_n,:);
    Back=p(B_n,:);
    Bottom=p(Bo_n,:);
    Top=p(T_n,:);
    %Right=p(L_n,:);
    %Left=p(R_n,:);
    FF=zeros(m*o,m*o);
    FB=DEA_SS(Front,Back,m*o,m*o,Beta);
    FBo=DEA_SS(Front,Bottom,m*o,n*m,Beta);
    FT=DEA_SS(Front,Top,m*o,n*m,Beta);
   
    BF=DEA_SS(Back,Front,m*o,m*o,Beta);
    BB=zeros(m*o,m*o);
    BBo=DEA_SS(Back,Bottom,m*o,n*m,Beta);
    BT=DEA_SS(Back,Top,m*o,n*m,Beta);

    BoF=DEA_SS(Bottom,Front,n*m,m*o,Beta);
    BoB=DEA_SS(Bottom,Back,n*m,m*o,Beta);
    BoBo=zeros(n*m,n*m);
    BoT=DEA_SS(Bottom,Top,n*m,n*m,Beta);
    
    TF=DEA_SS(Top,Front,n*m,m*o,Beta);
    TB=DEA_SS(Top,Back,n*m,m*o,Beta);
    TBo=DEA_SS(Top,Bottom,n*m,n*m,Beta);
    TT=zeros(n*m,n*m);
    
   
   

    SS=[FF FB FBo FT; BF BB BBo BT; BoF BoB BoBo BoT ; TF TB TBo TT];

function [SS_DEA_SS2]=DEA_SS(AA, BB, ii, jj,Beta)
SS_DEA_SS2=zeros(ii,jj);    
node = 0;
    for p2=1:jj
        for p1=1:ii
        node = node + 1;
        R= sqrt((AA(p1,1)-BB(p2,1)).^2 + (AA(p1,2)-BB(p2,2)).^2 + (AA(p1,3)-BB(p2,3)).^2);
        cos1=(AA(p1,4)*sqrt((AA(p1,1)-BB(p2,1)).^2)+AA(p1,5)*sqrt((AA(p1,2)-BB(p2,2)).^2)+AA(p1,6)*sqrt((AA(p1,3)-BB(p2,3)).^2))./R;
        cos2=(BB(p2,4)*sqrt((AA(p1,1)-BB(p2,1)).^2)+BB(p2,5)*sqrt((AA(p1,2)-BB(p2,2)).^2)+BB(p2,6)*sqrt((AA(p1,3)-BB(p2,3)).^2))./R;
        SS_DEA_SS2(p1,p2)=((exp(-(Beta).*R)).*((cos1.*cos2.*(AA(p1,7)).*(BB(p2,7)))))./((R.^2).*(pi));
        end
    end
end

  end