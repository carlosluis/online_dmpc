function noise = rnd_noise(std_pos, std_vel)

noise = [normrnd(0, std_pos, 3, 1); normrnd(0, std_vel, 3, 1)];
