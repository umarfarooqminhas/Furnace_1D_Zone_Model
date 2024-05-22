function [Temp_stock,Temp_stock_mid,Temp_s,Temp_mid] = Update_Stock_Pos(Bo_n,Temp_s,Temp_mid,n_s,Temp_s_array,Temp_mid_array,w,pos,Temp_stock,Temp_stock_mid)
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
         Temp_stock(w)=Temp_s(pos);                                                      % Stock temperature array along the length of furnace for final plotting
    Temp_stock_mid(w)=Temp_mid(pos); 
end
