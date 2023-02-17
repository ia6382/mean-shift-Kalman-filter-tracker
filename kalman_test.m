N = 20; %40
q = 1;
r = 1;

%krivulja1
%v = linspace(5*pi, 0, N);
%x = cos(v).*v;
%y = sin(v).*v;

%krivulja2
%v = linspace(5*pi, 0, N);
%x = v;
%y = sin(v);

%krivulja3
v = linspace(-3, 3, N);
x = v;
y = normpdf(x,0,0.1);

    %postopek za modele:
    %syms T q;
    %NCV
    %F = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    %L = [0 0; 0 0; 1 0; 0 1]
    %RW
    %F = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    %L = [1 0; 0 1; 0 0; 0 0]
    %NCA
    %F = [0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 0];
    %L = [0 0; 0 0; 0 0; 0 0; 1 0; 0 1];
    %A = expm(F*T); %Fi
    %Q=int((A*L)*q*(A*L)',T,0,T);

%ze dolocene (izpeljane) matrike:
%NCV
ncvA = [1 0 1 0;0 1 0 1; 0 0 1 0; 0 0 0 1]; %T = 1
ncvQ = q*[(1/3) 0 (1/2) 0; 0 (1/3) 0 (1/2); (1/2) 0 1 0; 0 (1/2) 0 1];
%RW
rwA = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
rwQ = q*[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0];
%NCA
ncaA = [1 0 1 0 (1/2) 0; 0 1 0 1 0 (1/2); 0 0 1 0 1 0; 0 0 0 1 0 1; 0 0 0 0 1 0; 0 0 0 0 0 1];
ncaQ = q*[(1/20) 0 (1/8) 0 (1/6) 0; 0 (1/20) 0 (1/8) 0 (1/6); (1/8) 0 (1/3) 0 (1/2) 0;  0 (1/8) 0 (1/3) 0 (1/2); (1/6) 0 (1/2) 0 1 0; 0 (1/6) 0 (1/2) 0 1];
ncaC = [1 0 0 0 0 0; 0 1 0 0 0 0];

C = [1 0 0 0; 0 1 0 0]; % za xyx.y.
R = r*[1 0; 0 1];

figure(1); 
subplot(1,3,1); kalman_plot(x, y, ncvA, ncvQ, C, R); pbaspect([1 1 1]); title(strcat('NCV: q= ',num2str(q),', r= ',num2str(r)));
%xlim([-1 17]);ylim([-1.2 1.2]);
subplot(1,3,2); kalman_plot(x, y, rwA, rwQ, C, R); pbaspect([1 1 1]); title(strcat('RW: q= ',num2str(q),', r= ',num2str(r)));
%xlim([-1 17]);ylim([-1.2 1.2]);
subplot(1,3,3); kalman_plot(x, y, ncaA, ncaQ, ncaC, R); pbaspect([1 1 1]); title(strcat('NCA: q= ',num2str(q),', r= ',num2str(r)));
%xlim([-1 17]);ylim([-1.2 1.2]);
hold off;