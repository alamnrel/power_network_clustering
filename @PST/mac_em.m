function [f] = mac_em(i,k,bus,flag)
% Syntax: [f] = mac_em(i,k,bus,flag)
%
% Purpose: generator electromechanical model, with vectorized computation
%          option. 
%          state variables are: mac_ang, mac_spd
%
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation and state state matrix
%                   building
%
% Output: f - dummy variable
%
% Files:
%
% See Also:

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     January 1991

% system variables
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot mac_int
global  mac_ang mac_spd eqprime edprime
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime

f = 0;

%% Initialization
jay = sqrt(-1);
if flag == 0;   
  [ngen, ncol] = size(mac_con);
  if ncol <= 21
    mac_con = [mac_con ones(ngen,23-ncol)];
  end
  if i ~= 0
    busnum = bus_int(mac_con(i,2)); % bus number
    mac_pot(i,1) = basmva/mac_con(i,3); % scaled MVA base
    mac_pot(i,2) = 1.0; % base kv
    % extract bus information
    eterm(i,1) = bus(busnum,2);   % terminal bus voltage
    theta(i,1) = bus(busnum,3)*pi/180;
    % terminal bus angle in radians
    pelect(i,1) = bus(busnum,4)*mac_con(i,22);
    % electrical output power, active
    qelect(i,1) = bus(busnum,5)*mac_con(i,23);
    % electrical output power, reactive
    curr = sqrt(pelect(i,1)^2+qelect(i,1)^2) ...
      /eterm(i,1)*mac_pot(i,1);  % current magnitude
    phi = atan2(qelect(i,1),pelect(i,1));
    % power factor angle
    v = eterm(i,1)*(cos(theta(i,1))+jay*sin(theta(i,1)));
    % voltage in real and imaginary parts
    % on system reference frame
    curr = curr*(cos(theta(i,1)-phi)+...
      jay*sin(theta(i,1)-phi)); % current in real and
    % imaginary parts on system reference frame
    eprime = v + jay*mac_con(i,7)*curr;
    ei = eprime;
    mac_ang(i,1) = atan2(imag(ei),real(ei));
    % machine angle (delta)
    mac_spd(i,1) = 1; % machine speed at steady state
    rot = sin(mac_ang(i,1))+jay*cos(mac_ang(i,1));
    % system reference frame rotation
    eprime = eprime*rot;
    edprime(i,1) = real(eprime);
    eqprime(i,1) = imag(eprime);
    curr = curr*rot;
    curdg(i,1) = real(curr); curqg(i,1) = imag(curr);
    curd(i,1) = real(curr)/mac_pot(i,1);
    curq(i,1) = imag(curr)/mac_pot(i,1);
    v = v*rot;
    ed(i,1) = real(v); eq(i,1) = imag(v);
    vex(i,1) = eqprime(i,1);
    pmech(i,1) = pelect(i,1)*mac_pot(i,1); % set input
    % mechanical power equal to electrical output power
    %keyboard
  else
    % vectorized computation
    [nmach,dum] = size(mac_con);
    busnum = bus_int(mac_con(:,2));  % bus number
    mac_pot1 = basmva*ones(nmach,1)./mac_con(:,3);  % scaled MVA base  % JHC 2013 1224
    mac_pot2 = 1.0*ones(nmach,1);  % base kv  % JHC 2013 1224
    mac_pot = [mac_pot1 mac_pot2];  % JHC 2013 1224
    % extract bus information
    eterm(:,1) = bus(busnum,2);         % terminal bus voltage
    theta(:,1) = bus(busnum,3)*pi/180;  % terminal bus angle in radians
    pelect(:,1) = bus(busnum,4).*mac_con(:,22);  % electrical output power, active
    qelect(:,1) = bus(busnum,5).*mac_con(:,23);  % electrical output power, reactive
    curr = sqrt(pelect(:,1).^2+qelect(:,1).^2) ...
      ./eterm(:,1).*mac_pot(:,1);  % current magnitude
    phi = atan2(qelect(:,1),pelect(:,1));  % power factor angle
    v = eterm(:,1).*(cos(theta(:,1))+jay*sin(theta(:,1)));
    % voltage in real and imaginary parts
    % on system reference frame
    curr = curr.*(cos(theta(:,1)-phi)+...
      jay*sin(theta(:,1)-phi)); % current in real and imaginary
    % parts on system reference frame
    eprime = v + jay*mac_con(:,7).*curr;
    ei = eprime;
    mac_ang(:,1) = atan2(imag(ei),real(ei));  % machine angle (delta)
    mac_spd(:,1) = ones(nmach,1);  % machine speed at steady state
    rot = sin(mac_ang(:,1))+jay*cos(mac_ang(:,1));  % system reference frame rotation
    eprime = eprime.*rot;
    edprime(:,1) = real(eprime);
    eqprime(:,1) = imag(eprime);
    curr = curr.*rot;
    curdg(:,1) = real(curr); curqg(:,1) = imag(curr);
    curd(:,1) = real(curr)./mac_pot(:,1);
    curq(:,1) = imag(curr)./mac_pot(:,1);
    v = v.*rot;
    ed(:,1) = real(v); eq(:,1) = imag(v);
    vex(:,1) = eqprime(:,1);
    pmech(:,1) = pelect(:,1).*mac_pot(:,1); % set input
    % mechanical power equal to electrical output power
  end
  %  keyboard
end

%% Network interface computation
if flag == 1 
  if i ~= 0
    mac_ang(i,k) = mac_ang(i,k) - mach_ref(k);   % wrt machine reference
    psi_re(i,k) = sin(mac_ang(i,k))*edprime(i,k) + ...
      cos(mac_ang(i,k))*eqprime(i,k);   % real part of psi
    psi_im(i,k) = -cos(mac_ang(i,k))*edprime(i,k) + ...
      sin(mac_ang(i,k))*eqprime(i,k);   % imag part of psi
    %keyboard
  else
    % vectorized computation
    [nmach,dum] = size(mac_con);
    mac_ang(:,k) = mac_ang(:,k)-mach_ref(k)*ones(nmach,1);   % wrt machine reference
    psi_re(:,k) = sin(mac_ang(:,k)).*edprime(:,k) + ...
      cos(mac_ang(:,k)).*eqprime(:,k);   % real part of psi
    psi_im(:,k) = -cos(mac_ang(:,k)).*edprime(:,k) + ...
      sin(mac_ang(:,k)).*eqprime(:,k);   % imag part of psi
  end
  %  keyboard
end

%% Generator dynamics calculation
if flag == 2 
  if i ~= 0
    curd(i,k) = sin(mac_ang(i,k))*cur_re(i,k) - ...
      cos(mac_ang(i,k))*cur_im(i,k); % d-axis current
    curq(i,k) = cos(mac_ang(i,k))*cur_re(i,k) + ...
      sin(mac_ang(i,k))*cur_im(i,k); % q-axis current
    curdg(i,k) = curd(i,k)*mac_pot(i,1);
    curqg(i,k) = curq(i,k)*mac_pot(i,1);
    dedprime(i,k) = 0;
    deqprime(i,k) = 0;
    ed(i,k) = edprime(i,k) + mac_con(i,7)*curqg(i,k);
    eq(i,k) = eqprime(i,k) - mac_con(i,7)*curdg(i,k);
    eterm(i,k) = sqrt(ed(i,k)^2+eq(i,k)^2);
    pelect(i,k) = eq(i,k)*curq(i,k) + ed(i,k)*curd(i,k);
    qelect(i,k) = eq(i,k)*curd(i,k) - ed(i,k)*curq(i,k);
    dmac_ang(i,k) = basrad*(mac_spd(i,k)-1.);
    dmac_spd(i,k) = (pmech(i,k)/mac_spd(i,k)-pelect(i,k)...
      *mac_pot(i,1)...
      -mac_con(i,17)*(mac_spd(i,k)-1))/(2.*mac_con(i,16));
    %keyboard
  else
    % vectorized computation
    [nmach,dum] = size(mac_con);
    curd(:,k) = sin(mac_ang(:,k)).*cur_re(:,k) - ...
      cos(mac_ang(:,k)).*cur_im(:,k); % d-axis current
    curq(:,k) = cos(mac_ang(:,k)).*cur_re(:,k) + ...
      sin(mac_ang(:,k)).*cur_im(:,k); % q-axis current
    curdg(:,k) = curd(:,k).*mac_pot(:,1);
    curqg(:,k) = curq(:,k).*mac_pot(:,1);
    dedprime(:,k) = zeros(nmach,1);
    deqprime(:,k) = zeros(nmach,1);
    ed(:,k) = edprime(:,k) + mac_con(:,7).*curqg(:,k);
    eq(:,k) = eqprime(:,k) - mac_con(:,7).*curdg(:,k);
    eterm(:,k) = sqrt(ed(:,k).^2+eq(:,k).^2);
    pelect(:,k) = eq(:,k).*curq(:,k) + ed(:,k).*curd(:,k);
    qelect(:,k) = eq(:,k).*curd(:,k) - ed(:,k).*curq(:,k);
    dmac_ang(:,k) = basrad*(mac_spd(:,k)-ones(nmach,1));
    while(0)
      dmac_spd(:,k) =(pmech(:,k)./mac_spd(:,k)...
        -pelect(:,k).*mac_pot(:,1)...
        -mac_con(:,17).*(mac_spd(:,k)...
        -ones(nmach,1)))./(2*mac_con(:,16));
    end %while(0)
    dmac_spd(:,k) =(pmech(:,k) ...
      -pelect(:,k).*mac_pot(:,1)...
      -mac_con(:,17).*(mac_spd(:,k)...
      -ones(nmach,1)))./(2*mac_con(:,16));
    
  end
  %  keyboard
end
