%Volume-Volume DEAs
function VV = get_volume_DEAs(vp,m,n,o,Beta)
    node = 0;
    VV=zeros(n*m*o,n*m*o);
    for vp2=1:(n*m*o)
        for vp1=1:(n*m*o)
        node = node + 1;
            if ((vp2~=vp1)) 
            R= sqrt((vp(vp1,1)-vp(vp2,1)).^2 + (vp(vp1,2)-vp(vp2,2)).^2 + (vp(vp1,3)-vp(vp2,3)).^2);
            VV(vp2,vp1)=(exp(-(Beta).*R).*(Beta.^2).*((vp(vp1,4)*vp(vp2,4))))./((R.^2).*(pi));
            else
            VV(vp2,vp1)=0;
            end
        end
    end
end



