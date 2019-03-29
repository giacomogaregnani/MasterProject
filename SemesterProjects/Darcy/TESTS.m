% Run test cases for DEM and CEM and one example of the DARCY PROBLEM
% References to the section 

% Choose the desired test
% For one dimensional case, set Test = 1 (Section 2.3.1)
% For two dimensional case, set Test = 2 (Section 2.3.2)
% For adaptivity vs DEM vs CEM, set Test = 3 (Seciton 2.3.3)
% For Darcy, set Test = 4 (Section 4)
clc; clear; close all;

Test = 3;

switch Test
    
    case 1
        run('OneDCase/Main.m')
        
    case 2
        run('TwoDCase/Main.m')
        
    case 3
        run('Adaptive2D/Main.m')
        
    case 4
        run('Darcy/Main.m')
        
end
        
