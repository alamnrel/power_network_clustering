function [T, i_best, contig] = best_part(cntg, all_bus, weXp, outl_gen, dif)

if nargin>=4
  min_gen_outl = min(outl_gen);
  gen_ind = outl_gen==min_gen_outl;
  
  if any(cntg(gen_ind))
    cntg(~gen_ind) = false;
    weXp(~cntg) = Inf;
    [~, i_best] = min(weXp);
    contig = 1;
  else
    dif(~gen_ind) = Inf;
    min_dif = min(dif);
    dif_ind = dif==min_dif;    
    weXp(~dif_ind) = Inf;
    weXp(~gen_ind) = Inf;
    [~, i_best] = min(weXp);    
    contig = 0;
  end
  
else % uncontsrained hMetis finds minimum worst expansion
  
  if any(cntg)
    weXp(~cntg) = Inf;
    [~, i_best] = min(weXp);
    contig = 1;
  else  % if no suitable partitions, choose most compact among most contiguous ones...
    min_dif = min(dif);
    dif_ind = dif==min_dif;
    weXp(~dif_ind) = Inf;
    [~, i_best] = min(weXp);
    contig = 0;
  end
end

T = all_bus{i_best};
end
