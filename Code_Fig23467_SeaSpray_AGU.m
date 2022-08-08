linspecer =...
   [63 64 250;...    blue... 
    238 51 56;...    red
    206 192 38;...   dark yellow
    187 125 180;...  magenta
    0 176 208;...    cyan
    0 138 77;...     green
    37 37 37;...     black
    196 122 47;...   brown
    245 132 53;...   orange
   128 128 128 ]/255;


load('u10.mat');
load('spray_volume.mat');
load('laser_ave_3min.mat');

P=laser134_avg_3mim;
R=26.5;
I0=1274.6;
Gamma=log(P./I0)/(-2*R);
V2=(1500.-I0.*exp(-2*R*Gamma)).*10^(-6)/3.75;
W=u10;


Toffoli_wind=(20:5:60);
Toffoli_wind_square=(20:5:60).^2;
Toffoli_wind_inverse=1.0./(20:5:60);
Toffoli_wind_inversesqrt=sqrt(1.0./(20:5:60));
Toffoli_laser0=[1500,1480,1450,1465,1360,1400,1268,1268,1100];
 Toffoli_laser=[1500,1480,1450,1405,1360,1314,1268,1184,1100];
Toffoli_laser_inverse=1.0./Toffoli_laser;
Toffoli_laser_inverse0=1.0./Toffoli_laser0;
% Toffoli_spray=[2.22e-7 5.36e-7 1.25e-6 0.92e-5 3.8e-5 4.57e-5  5.78e-5  7.187e-5 0.985e-4];
% Toffoli_spray_inverse=1.0./[2.22e-7 5.36e-7 1.25e-6 0.92e-5 3.8e-5 4.57e-5  5.78e-5  7.187e-5 0.985e-4];
Toffoli_spray=[2.22e-7 5.44e-7 1.21e-6 0.928e-5 3.75e-5 4.38e-5  5.55e-5 0.74e-4 0.97e-4];
Toffoli_spray_inverse=1.0./Toffoli_spray;
Toffoli_spray_log=log(Toffoli_spray);
Toffoli_spray_loginv=1.0./log(Toffoli_spray);

T_laser=Toffoli_laser(1:5);
T_laser0=Toffoli_laser0(1:5);
T_laser_inverse=Toffoli_laser_inverse(1:5);
T_laser_log=log(T_laser);
T_laser_log0=log(T_laser0);
T_wind=Toffoli_wind(1:5);
T_wind_inverse=Toffoli_wind_inverse(1:5);
T_spray=Toffoli_spray(1:5);
T_spray_log=log(T_spray);
% Gamma_field0=
Gamma_lab0=log(T_laser/1500)/(-2*1.29);
Gamma_Field=log(laser134_avg_3mim/1247.6)/(-2*26.5);
Spray_Wang=exp(136*Gamma_Field-15.25);
Spray_log=log(Spray_Wang);


%%
 % Following we reproduce the results in Toffoli et al. (2011),
 % Ma et al. (2020), and Xu et al. (2021)
 
 Laser_Field=laser134_avg_3mim;   % Laser Intensity in Field observation 
 Laser_Field2=Laser_Field;
%  Laser_Field2>

 
 Gamma_MaLab  =log(Toffoli_laser/1500)/(-2*1.3);   % Ma et al. (2020), laboratory P to ¦Ã
 Gamma_MaField=log(Laser_Field/1247.6)/(-2*26.5);  % Ma et al. (2021), field P to ¦Ã
 Gamma_XuLab  =log(Toffoli_laser/1742)/(-2*1.29)/2; % Xu et al. (2020), laboratory P to ¦Ì
 Gamma_XuField=log(Laser_Field2/1400)/(-2*26.5);   % Xu et al. (2020), Field P to ¦Ì

 
 %%  wind speed vs ¦Ã/¦Ì in Ma's and Xu's method

 %%%%%%%%%%%%%--------------------Figure 2------------------ %%%%%%%%%%%%
 figure;   % Wind speed vs ¦Ã/¦Ì£º Laboratory
 subplot(211)
 plot(Toffoli_wind,Gamma_MaLab,'o',Toffoli_wind, Gamma_XuLab,'o');
 xlabel('Wind in laboratory experiment (m s^{-1})'); ylabel('¦Ã or ¦Ì (m^{-1})');
 set(gca,'ytick',(0: 0.05: 0.2))
 xlim([15 65]); ylim([0 0.225])
 legend('¦Ã (Ma et al. 2020)','¦Ì (Xu et al. 2021)');
%  title('Laborotary Experiment');
 text(20, 0.2,'(a)'); grid on;
%  figure;  % Wind speed vs ¦Ã/¦Ì£º Field
 subplot(212)
 plot(u10,Gamma_MaField,'o',u10, Gamma_XuField,'o');
 xlabel('Wind in field experiment (m s^{-1})'); ylabel('¦Ã or ¦Ì (m^{-1})');
 legend('¦Ã (Ma et al. 2020)','¦Ì (Xu et al. 2021)');
%  title('Field Experiment');
 set(gca,'ytick',(0: 0.01: 0.04))
 set(gca,'xtick',(0: 5: 25),'xticklabel',num2str((0:5:25)'));
 xlim([0 25]); ylim([0 0.045])
 text(2.5, 0.04,'(b)');
 grid on;
 
%%
  
 % Spray volume flux in Ma et al. (2020): Field and Laboratory
 Spray_MaField=-5*10^(-3)*Gamma_MaField.*Gamma_MaField+1.5*10^(-3)*Gamma_MaField;
 Spray_MaLab  =-5*10^(-3)*Gamma_MaLab.*Gamma_MaLab+1.5*10^(-3)*Gamma_MaLab;
 
  % Spray volume flux in Xu et al. (2021): Field and Laboratory
 Spray_XuField=2.8*10^(3)*Gamma_XuField.^6.13;
 Spray_XuLab  =2.8*10^(3)*Gamma_XuLab.^6.13;

 
 XuLab=linspace(0.0,max(Gamma_XuLab),100);
%  XuLab=linspace(0.028,0.06,50);
 S_XuLab  =2.8*10^(3)*XuLab.^6.13;
 XuField=XuLab; XuField(XuField>max(Gamma_XuField))=[]; XuField(XuField<min(Gamma_XuField))=[];
 S_XuField=2.8*10^(3)*XuField.^6.13;

 MaLab=linspace(-0.01,max(Gamma_MaLab),100);
 S_MaLab  =-5*10^(-3)*MaLab.*MaLab+1.5*10^(-3)*MaLab;
 MaField=MaLab; MaField(MaField>max(Gamma_MaField))=[]; MaField(MaField<min(Gamma_MaField))=[];
 S_MaField=-5*10^(-3)*MaField.*MaField+1.5*10^(-3)*MaField;

  

 %%%%   ----------------Figure 3-----------------------   %%%%%%%%%%%
 figure;   % Wind speed vs ¦Ã/¦Ì£º Laboratory
 semilogy(Gamma_MaLab,Toffoli_spray,'o',Gamma_XuLab,Toffoli_spray,'o',...
           MaLab,S_MaLab,               XuLab, S_XuLab, ...
           MaField,S_MaField,'b',       XuField,S_XuField,'b');
 xlabel('¦Ã or ¦Ì (m^{-1})'); 
 ylabel('Spray volume flux (m^3 m^{-2} s^{-1})');
 set(gca,'xtick',(0:0.02:0.12))
 xlim([-0.01 0.13]); 
 ylim([10^(-8) 10^(-3)])
 legend('spray against ¦Ã for laboratory data','spray against ¦Ì for laboratory data',...
        'regression of spray onto ¦Ã','regression of spray onto ¦Ì');

    
    
    

%%%%%%%%%%% ----------------Figure 4 ---------------

figure;  % Laser Intensity vs Spray volume fluxï¼šToffoli
 err=[110, 35, 74, 60 40];
 h=errorbar(T_laser0,T_spray,err,'horizontal','*','LineWidth',1,'Color',[0 0.447 0.741]);
 h.MarkerSize=8;
 set(get(h,'Parent'), 'YScale', 'log')
   hold on;
 semilogy(T_laser0,T_spray,'*', ...
         T_laser,T_spray,'ko',...   
         (1355:1505),exp(-50*log((1355:1505))+350.5));
%  annotation('arrow',[T_laser0(4),T_spray(4)],[T_laser(4), T_spray(4)]);  
 legend('Errbar for laboratory data','laboratory data (Toffoli et al., 2011)','calibrated laboratory data',...
     'Regression of spray onto laser internsity')


 ylim([10^(-7) 10^(-4)])
 xlim([1300 1650])
 xlabel('Laser Internsity (W m^{-2})');
ylabel('spray volume flux (m^3 m^{-2} s^{-1})');

text(1365,3e-7,'R^2=0.99');
text(1365,2e-7,'RMSE=0.19');

% %%%%%%%%%%% ----------------Figure 5 ---------------
% figure;  % Wind speed vs Spray volume flux?
%  semilogy(u10,Spray_Wang,'o',(3:25),exp(-30./(3:25)-10.5));
%   ylim([10^(-7) 10^(-4)])
%  xlabel('Wind speed (m s^{-1})');
% ylabel('spray volume flux (m^3 m^{-2} s^{-1})');
% text(17,3e-7,'R^2=0.90');
% text(17,2e-7,'RMSE=0.30');


%%%%%%%%%%% ----------------Figure 6 ---------------
figure;  % Wind speed vs laser intensity in feild experiment
plot(u10,Laser_Field,'o',(3:25),400*exp(8./(3:25))-245);
  ylim([150 1500])
 xlabel('Wind speed (m s^{-1})');
ylabel('Laser intensity (W m^{-2})');
text(6,300,'R^2=0.94');
text(6,200,'RMSE=56.55');


%%%%%%%%%%% ----------------Figure 7 ---------------
figure;  % Wind speed vs Spray volume flux
semilogy(Toffoli_wind,Toffoli_spray,'o',(18:62),exp(-200./(18:62)-5.88));
ylim([10^(-7) 10^(-4)])
xlabel('Wind speed (m s^{-1})');
ylabel('Spray volume flux (m^3 m^{-2} s^{-1})');
text(50,3e-7,'R^2=0.95');
text(50,2e-7,'RMSE=0.54');

 
 

 