function[profile_new,Flux]=sinking(profile,z_int,omega,dt)


profile_new = profile;
thickness=z_int(2:end)-z_int(1:end-1);

Flux=profile.*omega.*dt;

profile_new(:)=profile(:)-Flux./thickness;
profile_new(2:end)=profile_new(2:end)+Flux(1:end-1)./thickness(2:end);
