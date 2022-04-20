%% ZSS PT1
data = readmatrix('../build/examples/PT1-Test.csv');
ct = data(:, 1);
cu = data(:, 2);
cy = data(:, 3);

K = 2;
T = 3;
Ts = 0.01;

s = tf('s');
Fs = K / (T * s + 1);
Fz = c2d(Fs, Ts);

[my1, mt1] = step(Fs);
[my2, mt2] = step(Fz);

close all;
figure
hold on
grid on
plot(mt1,my1)
stairs(mt2, my2)
stairs(ct, cy);
legend("Contious", "Discrete", "C")
title("PT1")


%% ZSS HPF
data = readmatrix('../build/examples/HPF-Test.csv');
ct = data(:, 1);
cu = data(:, 2);
cy = data(:, 3);

T = 1;
Ts = 0.1;

s = tf('s');
Fs = T * s / (T * s + 1);

z = zpk('z', Ts);
Fz = (z - 1) / ( z+ Ts/T - 1);

[my1, mt1] = step(Fs);
[my2, mt2] = step(Fz);

close all;
figure
hold on
grid on
plot(mt1,my1)
stairs(mt2, my2)
stairs(ct, cy);
legend("Contious", "Discrete", "C")
title("HPF")



%% ZSS PID v1
data = readmatrix('../build/examples/PID_v1-Test.csv');
ct = data(:, 1);
cw = data(:, 2);
ce = data(:, 3);
cy = data(:, 4);
cx = data(:, 5);

Ts = 0.1;

% PT1
K = 2;
T = 3;
s = tf('s');
Fs = K / (T * s + 1);
Fz = ss(c2d(Fs, Ts));
Fz.InputName = {'y'};
Fz.OutputName = {'x'};

% PID
Fpid = pid(1, 1, 0, 0.1, Ts);
Fpid.InputName = {'e'};
Fpid.OutputName = {'y'};

% System
Fw = connect(sumblk('e = w - x'), Fpid, Fz, {'w'}, {'w', 'e', 'y', 'x'});
[mv, mt] = step(Fw);
mw = mv(:, 1);
me = mv(:, 2);
my = mv(:, 3);
mx = mv(:, 4);


close all;
figure

% plot w
subplot(4, 1, 1)
hold on
grid on
stairs(mt, mw);
stairs(ct, cw);
legend("Matlab", "C++")
title("w")

% plot e
subplot(4, 1, 2)
hold on
grid on
stairs(mt, me);
stairs(ct, ce);
legend("Matlab", "C++")
title("e")

% plot y
subplot(4, 1, 3)
hold on
grid on
stairs(mt, my);
stairs(ct, cy);
legend("Matlab", "C++")
title("y")

% plot x
subplot(4, 1, 4)
hold on
grid on
stairs(mt, mx);
stairs(ct, cx);
legend("Matlab", "C++")
title("x")

%% ZSS PID v2
data = readmatrix('../build/examples/PID-Test_v2.csv');
ct = data(:, 1);
cw = data(:, 2);
ce = data(:, 3);
cy = data(:, 4);
cx = data(:, 5);

Ts = 0.1;

% HPF
T = 2;
z = zpk('z', Ts);
Fz = (z - 1) / ( z+ Ts/T - 1);
Fz = ss(Fz);
Fz.InputName = {'w'};
Fz.OutputName = {'e'};

% PID
Fpid = pid(0, 1, 0, 0.1, Ts);
Fpid.InputName = {'e'};
Fpid.OutputName = {'y'};

% System
Fw = connect(Fpid, Fz, {'w'}, {'w', 'e', 'y'});
[mv, mt] = step(Fw);
mw = mv(:, 1);
me = mv(:, 2);
my = mv(:, 3);


close all;
figure

% plot w
subplot(3, 1, 1)
hold on
grid on
stairs(mt, mw);
stairs(ct, cw);
legend("Matlab", "C++")
title("w")

% plot e
subplot(3, 1, 2)
hold on
grid on
stairs(mt, me);
stairs(ct, ce);
legend("Matlab", "C++")
title("e")

% plot y
subplot(3, 1, 3)
hold on
grid on
stairs(mt, my);
stairs(ct, cy);
legend("Matlab", "C++")
title("y")
