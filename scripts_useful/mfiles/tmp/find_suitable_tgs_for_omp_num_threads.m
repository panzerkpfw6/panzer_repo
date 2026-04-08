clear all
close all

function r = remain(nt,t_dim)
    r=rem(nt-2,(t_dim+1)*2);
end
