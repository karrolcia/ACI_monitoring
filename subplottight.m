function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    ax = subplot('Position', [(c-0.9)/m, 1-(r)/n, 0.75/m, 0.6/n]);
    if(nargout > 0)
      h = ax;
    end