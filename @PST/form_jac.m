function [Jac11,Jac12,Jac21,Jac22] = ...
      form_jac(V,ang,Y,ang_red,volt_red)
% Syntax:  [Jac] = form_jac(V,ang,Y,ang_red,volt_red)
%          [Jac11,Jac12,Jac21,Jac22] = form_jac(V,ang,Y,...
%                                      ang_red,volt_red)
%
% Purpose: form the Jacobian matrix
%
% Input:   V        - magnitude of bus voltage
%          ang      - angle(rad) of bus voltage
%          Y        - admittance matrix
%          ang_red  - vector to eliminate swing bus entries
%          volt_red - vector to eliminate generator bus
%                       entries
% Output:  Jac      - jacobian matrix
%          Jac11,Jac12,Jac21,Jac22 - submatrices of 
%                                      jacobian matrix  
%                                       
% See also:   
%
% Calls:
%
% Call By:   loadflow

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   1.0
% Author:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%
% ***********************************************************
jay = sqrt(-1);
[k dum] = size(Y);
cosang = cos(ang'); sinang = sin(ang');
% voltage perturbation rectangular coordinates
V_pert = cosang+jay*sinang;
% Voltage rectangular coordinates
V_rect = V'.*V_pert;
% angle and voltage perturbation rectangular coordinates
ang_pert = -V'.*(sinang-jay*cosang);
V_1 = conj(Y*V_rect);
% V_2 = diag(V_rect)*conj(Y);
%for i = 1:k,
%  V_2(i,:) = V_rect(i)*conj(Y(i,:));
%end
% sparse matrix formulation of V_2
i = [1:1:k]';
temp = sparse(i,i,V_rect,k,k);
V_2 = temp*conj(Y);

% X_1 = diag(V_1.*ang_pert)+V_2*diag(conj(ang_pert));
%for i = 1:k,
%  temp(:,i) = V_2(:,i)*conj(ang_pert(i));
%end
%X_1 = diag(V_1.*ang_pert) + temp;
% sparse matrix formulation of X_1
X_1 = sparse(i,i,V_1.*ang_pert,k,k);
X_1 = X_1 + V_2*sparse(i,i,conj(ang_pert),k,k);

%for i = 1:length(ang_red)
%  XX_1(:,i) = X_1(:,round(ang_red(i)));
%end
% sparse matrix formulation of XX_1
lang = length(ang_red);
ilang = [1:1:lang]';
x_red = sparse(round(ang_red),ilang,ones(lang,1),k,lang);
XX_1 = X_1*x_red;

%clear temp;
% X_2 = diag(V_1.*V_pert)+V_2*diag(conj(V_pert));
%for i = 1:k,
%  temp(:,i) = V_2(:,i)*conj(V_pert(i));
%end
%X_2 = diag(V_1.*V_pert) + temp;
% sparse matrix formulation of X_2
X_2 = sparse(i,i,V_1.*V_pert,k,k);
X_2 = X_2 + V_2*sparse(i,i,conj(V_pert),k,k);

%for i = 1:length(volt_red)
%  XX_2(:,i) = X_2(:,round(volt_red(i)));
%end
% sparse matrix formulation of XX_2
lvolt = length(volt_red);
ilvolt = [1:1:lvolt]';
x_volt = sparse(round(volt_red),ilvolt,ones(lvolt,1),k,lvolt);
XX_2 = X_2*x_volt;

%clear V_pert ang_pert
%save mem Y
%clear Y V_rect

% form submatrices of the Jacobian matrix
%J11 = ang_red*real(X_1)*ang_red';
%J21 = volt_red*imag(X_1)*ang_red';
%J12 = ang_red*real(X_2)*volt_red';
%J22 = volt_red*imag(X_2)*volt_red';
%for i = 1:length(ang_red),
%  J11(i,:) = real(XX_1(round(ang_red(i)),:));
%  J12(i,:) = real(XX_2(round(ang_red(i)),:));
%end
%for i = 1:length(volt_red),
%  J21(i,:) = imag(XX_1(round(volt_red(i)),:));
%  J22(i,:) = imag(XX_2(round(volt_red(i)),:));
%end
% sparse matrix formulation of J
temp = sparse(ilang,round(ang_red),ones(lang,1),lang,k);
J11 = temp*real(XX_1);
J12 = temp*real(XX_2);
temp = sparse(ilvolt,round(volt_red),ones(lvolt,1),lvolt,k);
J21 = temp*imag(XX_1);
J22 = temp*imag(XX_2);
if nargout > 3
   Jac11 = J11; clear J11
   Jac12 = J12; clear J12
   Jac21 = J21; clear J21
   Jac22 = J22; clear J22
  else
   Jac11 = [J11 J12;
	       J21 J22];
end

