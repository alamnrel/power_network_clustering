function mpw = mp_del_oos(mpw)
% Matpower Delete Out-Of-Service (mp_del_oos):
% Expertly remove all out-of-service elements and isolated buses by using
% the MATPOWER's ext2int utility, and then converting back to the original
% indexing. 
% That is, the function acts as a wrapper around ext2int with a small extra 
% functionality. 

% MATPOWER's ext2int (to remove out-of-service elements in mpw etc.):
mpw = ext2int(mpw);
i2e = mpw.order.bus.i2e;  

% Preserve the orinal bus labels of buses that remain:
mpw.bus(:,1) = i2e(mpw.bus(:,1));
mpw.gen(:,1) = i2e(mpw.gen(:,1));
mpw.branch(:,1) = i2e(mpw.branch(:,1));
mpw.branch(:,2) = i2e(mpw.branch(:,2));
end