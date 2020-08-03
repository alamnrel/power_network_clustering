function mpw = set_gencosts_svc(mpw, idx_pqg, qlvl)

ng = size(mpw.gen,1);
if nargin<2
  idx_pqg = true(ng,1);
else
  idx_pqg = idx_pqg(:);
end
if nargin<3
  qlvl = 0.48;
end
idx_pqg_p = idx_pqg;
idx_pqg_q = idx_pqg;

gencost_p = mpw.gencost;
gencost_q = mpw.gencost;
gencost_p(idx_pqg_p,1) = 1;
gencost_q(idx_pqg_q,1) = 1;
gencost_p(idx_pqg_p,2:3) = 0;
gencost_q(idx_pqg_q,2:3) = 0;
gencost_p(idx_pqg_p,4) = 2;
gencost_q(idx_pqg_q,4) = 3;
gencost_p(idx_pqg_p,5:6) = 0;
gencost_q(idx_pqg_q,5:6) = 0;
gencost_p(idx_pqg_p,7) = qlvl*mpw.gen(idx_pqg_p,9);
gencost_q(idx_pqg_q,7) = qlvl*mpw.gen(idx_pqg_q,4);
gencost_p(idx_pqg_p,8) = 1000;
gencost_q(idx_pqg_q,8) = 1000;
gencost_q(idx_pqg_q,9) = mpw.gen(idx_pqg_q,4);
gencost_q(idx_pqg_q,10) = 10000;

if size(gencost_p,2)<size(gencost_q,2)
  gencost_p = [gencost_p, zeros(ng, size(gencost_q,2)-size(gencost_p,2))];
end
if size(gencost_q,2)<size(gencost_p,2)
  gencost_q = [gencost_q, zeros(ng, size(gencost_p,2)-size(gencost_q,2))];
end

mpw.gencost = [gencost_p; gencost_q];
end