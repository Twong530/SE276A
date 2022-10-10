function [N_node,N_element,DOF,Nodal_position,Constant,EFT,FE,BCS] = input_CA1()

% Define nodes
N_node = 6
N_element = 10; 
E = 3e7*ones(N_element,1); 
A = ones(N_element,1);
Constant=[E,A]; 
DOF = 2;

% Define geometries
L = 20;
H = 15;

% Define nodal position
Nodal_position = [0, 0;
                  L, 0;
                  2*L,0;
                  3*L,0;
                  L, H;
                  2*L,H;];

% Define element free table (element, node 1, node 2)

EFT = [1, 1, 2;
       2, 2, 3;
       3, 3, 4;
       4, 1, 5;
       5, 2, 5;
       6, 3, 5;
       7, 2, 6;
       8, 3, 6;
       9, 4, 6;
       10, 5, 6;];

% Describe external force
FE = zeros(N_node*2,1); 
FE(3,1) = 0;
FE(4,1) = 0;
FE(5,1) = 0;
FE(6,1) = 0;
FE(7,1) = 0;
FE(9,1) = 3;
FE(10,1) = -9;
FE(11,1) = 0;
FE(12,1) = -12;

% Boundary conditions
BCS = zeros(1,N_node*2);  
BCS(1, 1:2) = 1;
BCS(1, 8) = 1;

end



