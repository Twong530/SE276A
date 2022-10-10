% Initiation
clear variables
close all
clc
format long

% Input
[N_node,N_element,DOF,Nodal_position,Constant,EFT,FE,BCS] = input_CA1();

% Define Matrix (two degrees of freedom)
K_global = zeros(N_node*DOF,N_node*DOF);
F_global = zeros(N_node*DOF,1);
Axial_stress = zeros(N_element,1);

%% FEM Assembly procedure
for ele=1:N_element

    E = Constant(ele,1);
    A = Constant(ele,2); 

    % Use element freedom table
    node_1 = EFT(ele,2);
    node_2 = EFT(ele,3);
    
    % Define the local matrix for each element
    L = sqrt((Nodal_position(node_2, 1)-Nodal_position(node_1, 1))^2 + (Nodal_position(node_2, 2)-Nodal_position(node_1, 2))^2); %length of element
    
    % Local element stiffness matrix
    K_local = (A*E/L)*[1 0 -1 0; 
                       0 0 0 0; 
                       -1 0 1 0; 
                       0 0 0 0;]; 

    r = sqrt((Nodal_position(node_2, 1)-Nodal_position(node_1, 1))^2+(Nodal_position(node_2, 2)-Nodal_position(node_2, 1))^2);
    s = ((Nodal_position(node_2, 1)-Nodal_position(node_1, 1))/r);
    c = ((Nodal_position(node_2, 2)-Nodal_position(node_1, 2))/r);

    % Find rotation matrix
    L_transformation = [c,s,0,0; 
                        -s,c,0,0; 
                        0,0,c,s; 
                        0,0,-s,c;];
    
    % Find element-wise stiffness matrix in global coordinates
    K_ele = (L_transformation')*K_local*(L_transformation);
    
    % Assembly into global stiffness matrix
    K_global(2*node_1-1:2*node_1,2*node_1-1:2*node_1) = K_global(2*node_1-1:2*node_1,2*node_1-1:2*node_1) + K_ele(1:2,1:2);
    K_global(2*node_1-1:2*node_1,2*node_2-1:2*node_2) = K_global(2*node_1-1:2*node_1,2*node_1-1:2*node_1) + K_ele(1:2,3:4);
    K_global(2*node_2-1:2*node_2,2*node_1-1:2*node_1) = K_global(2*node_1-1:2*node_1,2*node_1-1:2*node_1) + K_ele(3:4,1:2);
    K_global(2*node_2-1:2*node_2,2*node_2-1:2*node_2) = K_global(2*node_1-1:2*node_1,2*node_1-1:2*node_1) + K_ele(3:4,3:4);
    
end
%%
% Apply BCS
[K_global_withBC,F_withBC] = BoundaryConditionS_CA1(N_node,K_global,FE,BCS);
%%
% Solve
u_global = inv(K_global_withBC) * F_withBC
%% 
Force_reaction = K_global*u_global

%% Find axial force
for ele = 1:N_element
    node_1 = EFT(ele,2);
    node_2 = EFT(ele,3);
    r = sqrt((Nodal_position(node_2, 1)-Nodal_position(node_1, 1))^2+(Nodal_position(node_2, 2)-Nodal_position(node_2, 1))^2);

    % Find the angle
    s = ((Nodal_position(node_2, 1)-Nodal_position(node_1, 1))/r);
    c = ((Nodal_position(node_2, 2)-Nodal_position(node_1, 2))/r);

    Axial_stress(ele,1) = Force_reaction(2*node_2-1,1)*c+Force_reaction(2*node_2,1)*s-Force_reaction(2*node_1-1,1)*c-Force_reaction(2*node_1,1)*s;
end

