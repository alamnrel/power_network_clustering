function mpw = struct_mp(bus, line, gen, basmva, gencost)

if nargin<5
  ng = size(gen,1);
  gencost = zeros(ng,7);
end

mpw = struct('version', '2', 'bus', bus, 'branch', line, 'gen', gen, 'baseMVA', basmva, 'gencost', gencost);
end