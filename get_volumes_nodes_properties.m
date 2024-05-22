function vp = get_volumes_nodes_properties(m,n,o,L,W,H)
node = 0;                                                                           % Volumes nodes information
Vol=zeros(m,n);
for k = 1:o
    for j = 1:m
        for i = 1:n
        node = node + 1;
        Vol(node,1) = W/n*(i-0.5);
        Vol(node,2) = L/m*(j-0.5);
        Vol(node,3) = H/o*(k-0.5);
        Vol(node,4) = (W/n)*(L/m)*(H/o);
        Vol(node,5) = (W)*(L)*(H);
        end
    end
end

vp=Vol;                                                                           % All Volume nodes information