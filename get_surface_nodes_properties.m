function p = get_surface_nodes_properties(m,n,o,L,W,H)
Front=S1(m,o,W,L,m,H,o);                                                            % Front wall Surface Nodes Data (x,y,z, normal vector, element area, entire surface area)
Back=S2(m,o,W,L,m,H,n,o);                                                           % Back walls nodes information
Bottom=S3(n,m,W,L,m,H,n,o);                                                         % Stock Nodes information
Top=S4(n,m,W,L,m,H,n);                                                              % Top surface Nodes information
%Left=S5(n,o,W,L,m,H,n,o);
%Right=S6(n,o,W,L,H,n,o);                                                           % Left and right surfaces considered as entry and exit point

p=[Front;Back;Bottom;Top];                                                          % All Surface nodes information


function Front=S1(jj, kk,W,L,m,H,o)
node = 0;
Front=zeros(jj,kk);
    for j = 1:jj
        for k = 1:kk
        node = node + 1;
        Front(node,1) = W; % Front Surface
        Front(node,2) = L/m*(j-0.5);
        Front(node,3) = H/o*(k-0.5);
        Front(node,4) = 1;
        Front(node,5) = 0;
        Front(node,7) = (L/m)*(H/o);
        Front(node,8) = (L)*(H);
        end
    end
end




function [Back]=S2(jj, kk,W,L,m,H,n,o)
Back=zeros(jj,kk);
node = 0;
    for j = 1:jj
        for k = 1:kk
        node = node + 1;
        Back(node,1) = W/n*0; % Back Surface
        Back(node,2) = L/m*(j-0.5);
        Back(node,3) = H/o*(k-0.5);
        Back(node,4) = 1;
        Back(node,5) = 0;
        Back(node,6) = 0;
        Back(node,7) = (L/m)*(H/o);
        Back(node,8) = (L)*(H);
        end
    end
end

function [Bottom]=S3(ii, jj,W,L,m,H,n,o)
node = 0;
Bottom=zeros(jj,ii);
    for j = 1:jj
        for i = 1:ii
        node = node + 1;
        Bottom(node,1) = W/n*(i-0.5); % Bottom surface
        Bottom(node,2) = L/m*(j-0.5);
        Bottom(node,3) = H/o*0;
        Bottom(node,4) = 0;
        Bottom(node,5) = 0;
        Bottom(node,6) = 1;
        Bottom(node,7) =(W/n)*(L/m);
        Bottom(node,8) = (W)*(L);
        end
    end
end

function [Top]=S4(ii, jj,W,L,m,H,n)
node = 0;
Top=zeros(jj,ii);
    for j = 1:jj
        for i = 1:ii
        node = node + 1;
        Top(node,1) = W/n*(i-0.5); % Top surface
        Top(node,2) = L/m*(j-0.5);
        Top(node,3) = H;
        Top(node,4) = 0;
        Top(node,5) = 0;
        Top(node,6) = 1;
        Top(node,7) = (W/n)*(L/m);
        Top(node,8) = (W)*(L);
        end
    end
end
%{
function [Left]=S5(ii, kk,W,L,m,H,n,o)
node = 0;
Left=zeros(ii,kk);
    for i = 1:ii
        for k = 1:kk
        node = node + 1;
        Left(node,1) = W/n*(i-0.5); % Left surface
        Left(node,2) = L/m*0;
        Left(node,3) = H/o*(k-0.5);
        Left(node,4) = 0;
        Left(node,5) = 1;
        Left(node,6) = 0;
        Left(node,7) = (W/n)*(H/o);
        Left(node,8) = (W)*(H);
        end
    end
end

function [Right]=S6(ii, kk,W,L,H,n,o)
node = 0;
Right=zeros(ii,kk);
    for i = 1:ii
        for k = 1:kk
        node = node + 1;
        Right(node,1) = W/n*(i-0.5); % Right surface
        Right(node,2) = L;
        Right(node,3) = H/o*(k-0.5);
        Right(node,4) = 0;
        Right(node,5) = 1;
        Right(node,6) = 0;
        Right(node,7) = (W/n)*(H/o);
        Right(node,8) = (W)*(H);
        end
    end
end
%}
end