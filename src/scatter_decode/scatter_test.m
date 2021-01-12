%% gen_IK
m = 4;
I = gen_IK(m)

%% gen_rk
r = exp(1j*2*pi/4*[1 2 1 1 1 2 0 1 3 0 2 3 3 0 2 1].');
rk = gen_rk(r);
mod(angle(rk), 2*pi) / (2*pi) * 4

%% gen_index

