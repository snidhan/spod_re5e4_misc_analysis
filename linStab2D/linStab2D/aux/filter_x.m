function [f] = filter_x(f,alpha_f)

[nx,ny]=size(f);

for i=1:ny
   f(2:end-1,i) =  filter_O10(f(2:end-1,i),alpha_f);
end

end

