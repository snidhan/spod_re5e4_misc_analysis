function [f] = filter_y(f,alpha_f)

[nx,ny]=size(f);

for i=1:nx
   f(i,2:end-1) =  filter_O10(f(i,2:end-1),alpha_f);
end

end

