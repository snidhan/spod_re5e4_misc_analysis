spod_mode  = spod_modes1_arranged(:,3,,1);
spod_mode_mag = spod_mode.*conj(spod_mode);
spod_mode_mag = spod_mode_mag.*rc;
eigvalue = trapz(rc,spod_mode_mag);
