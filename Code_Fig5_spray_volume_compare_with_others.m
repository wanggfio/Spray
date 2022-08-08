%##########################################################################################
%------------------------------------------------------------------------------------------
% --- Comparision with other works.
% --- field observation.
clc; clear; close all;

% linspecer =...
%    [63 64 250;...    blue... 
%     238 51 56;...    red
%     206 192 38;...   dark yellow
%     187 125 180;...  magenta
%     0 176 208;...    cyan
%     0 138 77;...     green
%     37 37 37;...     black
%     196 122 47;...   brown
%     245 132 53;...   orange
%    128 128 128 ]/255;

linspecer =...
   [0 0.447 0.741;...    
    0.85 0.325 0.098;...    
    0.929 0.694  0.125;...  
    0.494  0.184  0.556;...  
    0.466 0.674  0.188;...    
    0.301 0.745 0.933;...     
    0.635 0.078 0.184;...     
    0.769 0.478 0.184;...   
    0.961 0.518 0.208;...   
    0.145 0.145 0.145 ];

load('.\spray_volume.mat');

%------------------------------------------------------------------------------------------
% --- Andreas (1998).
load('.\andreas_1998.txt')
u10_andreas=andreas_1998(:,1);
SV_andreas=andreas_1998(:,2);
%------------------------------------------------------------------------------------------
% --- Toffoli et al. (2011).
load('.\Toffoli_et_al_2011_spray_volume.txt')
u10_toffoli=Toffoli_et_al_2011_spray_volume(:,1);
SV_toffoli=Toffoli_et_al_2011_spray_volume(:,2);
%------------------------------------------------------------------------------------------
% --- The work of Zhao et al.,(2006)
u10_zhao = 5 : 1 : 30;
beta = 0.2;
cd = 2*10^(-3);
v_air = 1.48*10^(-5);
rb = cd * beta .* u10_zhao.^3./(9.8 * v_air);
SV_zhao_low = 0;
for r = 30 : 1 : 75
    SV_zhao_low = SV_zhao_low + 7.84*10^(-3).*rb.^(1.5)*r^(-1)*4/3*3.14*(r*10^-6)^3;
end
for r = 76 : 1 : 200
    SV_zhao_low = SV_zhao_low + 4.41*10^(1).*rb.^(1.5)*r^(-3)*4/3*3.14*(r*10^-6)^3;
end
for r = 201 : 1 : 500
    SV_zhao_low = SV_zhao_low + 1.41*10^(13).*rb.^(1.5)*r^(-8)*4/3*3.14*(r*10^-6)^3;
end

u10_zhao = 5 : 1 : 30;
beta = 1.2;
cd = 2*10^(-3);
v_air = 1.48*10^(-5);
rb = cd * beta .* u10_zhao.^3./(9.8 * v_air);
SV_zhao_high = 0;
for r = 30 : 1 : 75
    SV_zhao_high = SV_zhao_high + 7.84*10^(-3).*rb.^(1.5)*r^(-1)*4/3*3.14*(r*10^-6)^3;
end
for r = 76 : 1 : 200
    SV_zhao_high = SV_zhao_high + 4.41*10^(1).*rb.^(1.5)*r^(-3)*4/3*3.14*(r*10^-6)^3;
end
for r = 201 : 1 : 500
    SV_zhao_high = SV_zhao_high + 1.41*10^(13).*rb.^(1.5)*r^(-8)*4/3*3.14*(r*10^-6)^3;
end
%------------------------------------------------------------------------------------------
% --- Fairall et al. (1994)
u10 = [6 9 11 13 15 18];
B0 = [ 3.726  4.138  4.405  4.596  4.758  4.998];
B1 = [-3.656 -3.236 -2.646 -2.232 -2.038 -1.758];
B2 = [ 3.673  1.172 -3.156 -5.983 -7.101 -9.323];
B3 = [-0.629  2.292  8.902  13.198 14.758 18.238];
B4 = [-0.525 -1.569 -4.482 -6.382 -7.038 -8.403];
C1 = [4.52*10^3 7.18*10^3 1.02*10^4 1.36*10^4 1.81*10^4 6.34*10^4];
C2 = [3.08*10^6 4.89*10^6 6.95*10^6 9.25*10^6 1.23*10^7 4.32*10^7];
C3 = [7.73*10^16 1.23*10^17 1.75*10^17 2.32*10^17 3.10*10^17 1.08*10^18];
df_num = 0;
for ii = 1 : 1 :6
    num = 1;
    for r80 = 9.6 : 0.1 : 15
        r = (r80/0.518)^(1/0.976);
        dfdr(num) = 0.506*r^(-0.024)*(10^(B0(ii) + B1(ii)*log10(r80) + B2(ii)*(log10(r80))^2 + B3(ii)*(log10(r80))^3 + B4(ii)*(log10(r80))^4 ));
        num = num +1;
    end
    for r80 = 15.1 : 0.1 : 37.5
        r = (r80/0.518)^(1/0.976);
        dfdr(num) = 0.506*r^(-0.024)*(C1(ii)*r80^(-1));
        num = num +1;
    end
    for r80 = 37.6 : 0.1 : 100
        r = (r80/0.518)^(1/0.976);
        dfdr(num) = 0.506*r^(-0.024)*(C2(ii)*r80^(-2.8));
        num = num +1;
    end
    for r80 = 100.1 : 0.1 : 200
        r = (r80/0.518)^(1/0.976);
        dfdr(num) = 0.506*r^(-0.024)*(C3(ii)*r80^(-8));
        num = num +1;
    end
    df_num = df_num + dfdr./(3.8*10^-6*u10(ii)^3.4);
end

% --- Fairall et al.(1994), fn(r).
    df_num_avg = df_num./6;


% --- zhao et al.(2006),the dash line of fairall et al.(1994) in figure 8.
wind = 10 : 1 : 30;
for ii = 1 : 21
    DF(ii) = sum(3.8 * 10^-6*(wind(ii))^3.4 .* df_num_avg.*0.1);
end



%figure;
%semilogy(wind, DF)
u10_Fairall = 5 : 1 : 30;
r80 = 9.6 : 0.1 : 200;
SV_Fairall(1 :  size(u10_Fairall,2)) = 0;
for jj = 1 : 1 : size(u10_Fairall,2)
    for ii = 1 : 1 : size(r80,2)
        r = (r80(ii)/0.518)^(1/0.976);
        SV_Fairall(jj) = SV_Fairall(jj) + 3.8 * 10^-6*(u10_Fairall(jj))^3.4 * df_num_avg(ii)*0.1*(4/3*3.14*(r*10^-6)^3);
    end
end
%------------------------------------------------------------------------------------------
% --- Monahan et al. 1986.
u10_monahan = 5 : 0.5 : 30 ;
SV_monahan(1 : size(u10_monahan,2))=0;
for jj = 1 : 1 : size(u10_monahan,2)
    for r = 1 : 1 : 10
        dfdr_monahan(r) = 1.373*u10_monahan(jj)^3.41*r^(-3)*(1+0.057*r^1.05)*10^(1.19*exp(-((0.380-log10(r))/0.650)^2))*(4/3*3.14*(r*10^-6)^3);
    end
    for r = 11 : 1 : 75
        dfdr_monahan(r) = 8.60*10^(-6)*exp(2.08*u10_monahan(jj))*r^(-2)*(4/3*3.14*(r*10^-6)^3);
    end
    for r = 76 : 1 : 100
        dfdr_monahan(r) = 4.83*10^(-2)*exp(2.08*u10_monahan(jj))*r^(-4)*(4/3*3.14*(r*10^-6)^3);
    end
    for r = 101 : 1 : 500
        dfdr_monahan(r) = 4.83*10^6*exp(2.08*u10_monahan(jj))*r^(-8)*(4/3*3.14*(r*10^-6)^3);
    end
    SV_monahan(jj) = sum(dfdr_monahan);
end

%------------------------------------------------------------------------------------------
% --- Iida et al.,(1992);
u10_Iida = 5 : 30;
cd = 2*10^(-3);
v_air = 1.48*10^(-5);
u_star_sq = cd.*u10_Iida.*u10_Iida;

c0 = [2.8 2.8 2.8 2.8 2.8 2.8 2.8 2.8];
c1 = [-12.25 -12.33 -12.11 -12.17 -12.16 -12.45 -12.48 -12.79];
c0_2 = [3.90 3.90 3.90 3.90 3.90 3.90 ];
c1_2 = [-17.22 -17.27 -16.81 -16.45 -16.43 -16.57];
dr = [13 21 29 42 61 88 130 191];

for j = 1 : 26
    u_star_sq = cd.*u10_Iida(j).*u10_Iida(j);

    beta = 0.2;
    if u_star_sq.*u10_Iida(j)/v_air/9.8*beta >= 10^4
        for i = 1 : 8
            theta(i) = 10^(c0(i) * log10(u_star_sq.*u10_Iida(j)/v_air/9.8*beta) + c1(i))/dr(i)*10^3;
            %theta(i) = 10^(c0(i) * log10(2.24*10^4) + c1(i))/dr(i)*10^3;
        end
        radius = 28 : 595;
        df(1 : 13) = theta(1);df(14 : 13+20) = theta(2);df(13+20+1 : 61) = theta(3);df(62 : 102) = theta(4);
        df(103 : 162) = theta(5);df(163 : 249) = theta(6);df(250 : 378) = theta(7);df(379 : 568) = theta(8);
        df_beta_low = df;
        sum_low(j) =0;
        for i = 1 : 568
            sum_low(j) = sum_low(j) + df_beta_low(i)*(4/3*pi*(radius(i)*10^-6 )^3);
        end
    else
        for i = 1 : 6
            theta(i) = 10^(c0_2(i) * log10(u_star_sq.*u10_Iida(j)/v_air/9.8*beta) + c1_2(i))/dr(i)*10^3;
            %theta(i) = 10^(c0(i) * log10(2.24*10^4) + c1(i))/dr(i)*10^3;
        end
        radius = 28 : 276;
        df(1 : 13) = theta(1);df(14 : 13+20) = theta(2);df(13+20+1 : 61) = theta(3);df(62 : 102) = theta(4);
        df(103 : 162) = theta(5);df(163 : 249) = theta(6);
        df_beta_low = df;
        sum_low(j) =0;
        for i = 1 : 249
            sum_low(j) = sum_low(j) + df_beta_low(i)*(4/3*pi*(radius(i)*10^-6 )^3);
        end        
    end

    beta = 1.2;
    if u_star_sq.*u10_Iida(j)/v_air/9.8*beta >= 10^4
        for i = 1 : 8
            theta(i) = 10^(c0(i) * log10(u_star_sq.*u10_Iida(j)/v_air/9.8*beta) + c1(i))/dr(i)*10^3;
            %theta(i) = 10^(c0(i) * log10(1.34*10^5) + c1(i))/dr(i)*10^3;
        end
        radius = 28 : 595;
        df(1 : 13) = theta(1);df(14 : 13+20) = theta(2);df(13+20+1 : 61) = theta(3);df(62 : 102) = theta(4);
        df(103 : 162) = theta(5);df(163 : 249) = theta(6);df(250 : 378) = theta(7);df(379 : 568) = theta(8);
        df_beta_high = df;
        sum_high(j) =0;
        for i = 1 : 568
            sum_high(j) = sum_high(j) + df_beta_high(i)*(4/3*pi*(radius(i)*10^-6 )^3);
        end
    else
        for i = 1 : 6
            theta(i) = 10^(c0_2(i) * log10(u_star_sq.*u10_Iida(j)/v_air/9.8*beta) + c1_2(i))/dr(i)*10^3;
            %theta(i) = 10^(c0(i) * log10(1.34*10^5) + c1(i))/dr(i)*10^3;
        end
        radius = 28 : 276;
        df(1 : 13) = theta(1);df(14 : 13+20) = theta(2);df(13+20+1 : 61) = theta(3);df(62 : 102) = theta(4);
        df(103 : 162) = theta(5);df(163 : 249) = theta(6);
        df_beta_high = df;
        sum_high(j) =0;
        for i = 1 : 249
            sum_high(j) = sum_high(j) + df_beta_high(i)*(4/3*pi*(radius(i)*10^-6 )^3);
        end        
    end
end
%figure(1);
%loglog(radius, df_beta_low);hold on;
%loglog(radius, df_beta_high);hold off;
%ylim([10^-3 10^5])
%grid on;
%------------------------------------------------------------------------------------------


% --- Troitskaya et al. 2018 Part I.
load('.\spray_vlume_troitskaya_2018a.mat')

%%
load('.\u10.mat');
load('.\laser_ave_3min.mat');
Gamma_Field=log(laser134_avg_3mim/1247.6)/(-2*26.5);
Spray_Wang=exp(136*Gamma_Field-15.25);
%  Laser_Field=laser134_avg_3mim;
%  Gamma_XuField=log(Laser_Field/1400)/(-2*26.5); 
%  Spray_XuField=2.8*10^(3)*Gamma_XuField.^6.13;
figure(8);
set(gcf,'position',[200 100 600 800]);

semilogy(u10_monahan,SV_monahan,'kd');hold on;
semilogy(u10_Iida(2 : end), sum_low(2 : end),'^','MarkerEdgeColor',linspecer(4,:));hold on;
semilogy(u10_Iida(2 : end), sum_high(2 : end),'^','MarkerEdgeColor',linspecer(5,:));hold on;
semilogy(u10_Fairall(2 : end), SV_Fairall(2 : end),'v','MarkerEdgeColor',linspecer(6,:));hold on;
semilogy(u10_andreas, SV_andreas, 'o','MarkerEdgeColor',linspecer(7,:));hold on;
semilogy(u10_zhao(6 : end), SV_zhao_low(6 : end),'s','MarkerEdgeColor',linspecer(8,:));hold on;
semilogy(u10_zhao(6 : end), SV_zhao_high(6 : end),'s','MarkerEdgeColor',linspecer(9,:));hold on;
semilogy(u10_toffoli, SV_toffoli,'.','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',20);hold on;
semilogy(u10_troi, SV_troi,'.','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',20);hold on;
% semilogy(u10,Spray_Wang,'o');hold on;
semilogy(u10, V_NR,'k.');hold on;       % Ma et al., 2019
% semilogy((5:20), -14*linspace(exp(5),exp(20),length(5:20)),'ro');    % Xu et al., 2021
semilogy([5 20], (1.8e-8)*[5 20]-(6e-8),'-k'); hold on;   % Xu et al., 2021
semilogy(u10,Spray_Wang,'o','MarkerEdgeColor',linspecer(1,:));hold on;
semilogy((4:30),exp(-30./(4:30)-10.5),'color',linspecer(2,:));
%semilogy(windspeed_14p77_avg_3mim,spray_volume_cubic,'r.');hold off;
%legend('Considering Rangefinder equation', 'No considering Rangefinder equation','Location','Southeast')
xlim([0 36]);
ylim([10^-11 10^-3 ])
grid on;
%title('Wind speed at 14.77m height')
legend1 = legend('Monahan et al. (1986)','Iida et al. (1992): \beta = 0.2','Iida et al. (1992): \beta = 1.2', ... 
                 'Fairall et al. (1994)','Andreas (1998)', 'Zhao et al. (2006): \beta = 0.2','Zhao et al. (2006): \beta = 1.2',...
                 'Toffoli et al. (2011)','Troitskaya et al. (2018)','Ma et al. (2020)','Xu et al. (2021)','This study','This study: Regression', 'Location','Southeast');
set(legend1,'fontweight','bold','Fontname', 'Times New Roman');
set(legend1,'fontsize',12);
% legend('boxoff')
xlabel('Wind speed (m s^-^1)','fontsize',14,'fontweight','bold','Fontname', 'Times New Roman');
ylabel('Spray volume flux (m^3 m^-^2 s^-^1)','fontsize',14,'fontweight','bold','Fontname', 'Times New Roman')
set(gca,'fontsize',12);set(gca,'fontsize',14,'fontweight','bold','Fontname', 'Times New Roman');
set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 4.5 6]);
print(gcf,'-dpng','-r1500','spray_volume_comp_1');
%%