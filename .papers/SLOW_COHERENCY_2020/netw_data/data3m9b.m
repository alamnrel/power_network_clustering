% A 3-machine 9-bus system from Chow's book pp.70
% data3m9b.m

% bus data format
% bus: number, voltage(pu), angle(degree), p_gen(pu), q_gen(pu),
%      p_load(pu), q_load(pu), bus_type
%      bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)

bus = [ 1 1.04    0.00   0.00   0.00  0.00  0.00  0.00  0.00 1;
	2 1.02533 0.00   1.63   0.00  0.00  0.00  0.00  0.00 2;
	3 1.02536 0.00   0.85   0.00  0.00  0.00  0.00  0.00 2;
	4 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00 3;
	5 1.00    0.00   0.00   0.00  0.90  0.30  0.00  0.00 3;
	6 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00 3;
	7 1.00    0.00   0.00   0.00  1.00  0.35  0.00  0.00 3;
	8 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00 3;
	9 1.00    0.00   0.00   0.00  1.25  0.50  0.00  0.00 3];

% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio

line = [1 4 0.0    0.0576 0.     1. 0. ;
	4 5 0.017  0.092  0.158  1. 0. ;
	5 6 0.039  0.17   0.358  1. 0. ;
	3 6 0.0    0.0586 0.     1. 0. ;
	6 7 0.0119 0.1008 0.209  1. 0. ;
	7 8 0.0085 0.072  0.149  1. 0. ;
	8 2 0.0    0.0625 0.     1. 0. ;
	8 9 0.032  0.161  0.306  1. 0. ;
	9 4 0.01   0.085  0.176  1. 0. ];

% Machine data (The data is from Anderson-Fouad! Thus the inertias are 
%  lower than in Chow-1982, but the rest is the same)
mac_con = [1,2,3;
           1,2,3;
           100,100,100;
           0.0336,0.0521,0.0742;
           0,0,0;
           0.146,0.896,1.312;
           0.0608,0.12,0.18;
           0,0,0;
           8.96,6.00,5.89;
           0,0,0;
           0.097,0.864,1.258;
           0.097,0.197,0.25;           
           0,0,0;
           0,0.535,0.6;
           0,0,0;
           %2364/247.5, 640/192, 301/128;   %inertias as in Anderson
           2364/100, 640/100, 301/100;   %inertias as in Chow
           0,0,0;
           0,0,0;
           1,2,3];         
mac_con = mac_con';