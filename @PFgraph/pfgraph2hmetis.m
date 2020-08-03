function pfgraph2hmetis(g, vw)
%
% The node weights need to be provided separately, as PFgraph node weights
% are meant for optimizing the normalized cut problem (g.adj, g.vw). On the
% other hand, the goal of using HMETIS may be entirely different. If HMETIS
% is used with node weights for normalized spectral clustering, this needs
% to be re-confirmed by setting the vw input as g.vw.
%

if nargin>1
  GraphUtils.pfgraph2hmetis(g.adj, vw);
else
  GraphUtils.pfgraph2hmetis(g.adj);
end


end