y0 = [1; 2];
t_span = [0, 100];
[t, y] = ode45(@example_ode1, t_span, y0);