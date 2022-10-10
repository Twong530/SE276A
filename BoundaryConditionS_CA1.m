function [K_global_withBC,F_withBC] = BoundaryConditionS_CA1(N_node,K_global,FE,BCS)
%   Define boundary conditions and get K reduced  
K_global_withBC = K_global;
F_withBC = FE;

for i = 1:2*N_node
        if BCS(1,i) == 1
             K_global_withBC(i,:) = 0;
             K_global_withBC(:,i) = 0;
             K_global_withBC(i,i) = 1;
             F_withBC(i,1) = 0; 
             %only for case with 0 prescribed BC
        end
end
end