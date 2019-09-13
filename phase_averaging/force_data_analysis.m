%%

clear;
filename = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/force_files/forces_xyz_imb0001_t913_t1374.plt';
A = importdata(filename);
time = A(:,1);

time_final(1,1) = time(1,1);
count = 1;

for i = 2:size(time,1)-1
    t1 =  time_final(count,1);
    t2 =  time(i+1,1);
    if t2 > t1
        count = count + 1;
        time_final(count,1)  = time(i+1,1);
        A_final(count,:) = A(i+1,:);
    end 
end


%%

fz = A_final(223:end,end);

time_start = time_final(223,1);
time_end = time_final(end,1);


size_time_signal = size(time_final,1) - 223 + 1;
time_uniform = linspace(time_start, time_end, size_time_signal)';

fz_uniform = spline(time_final(223:end,1), fz, time_uniform);
fz_uniform = fz_uniform';
fz_uniform_fluc = fz_uniform - mean(fz_uniform);
%% 

hilb_fz_uniform_fluc = hilbert(fz_uniform_fluc);
analytic_fz = fz_uniform_fluc - sqrt(-1)*hilb_fz_uniform_fluc;
phase_angle = atan(imag(analytic_fz)./real(analytic_fz));

plot(time_uniform,phase_angle,'ko');