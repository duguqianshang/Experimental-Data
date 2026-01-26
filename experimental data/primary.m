clc
clear
close all
prgb=readmatrix('');
prgs=readmatrix('');
prhb=readmatrix('');
PV1=readmatrix('');
WT1=readmatrix('');
PL10=readmatrix('');
HL10=readmatrix('');
GL10=readmatrix('');
LB1_recoder=0;
T=192;
Png_chp1_max=1000;
Hng_bo1_max=700;
Esp1_max=50;
GP_max=1000;
GH_max=600;
eff_ch=0.95;
eff_dch=0.95;
D1=0.8;
M=10000;
K=548;
pv_ybsj=readmatrix(''); 
wt_ybsj1=readmatrix(''); 
wc=(PV1+WT1)-(pv_ybsj+wt_ybsj1)';
wc=wc';
[wc_xj,wc_sj,sj_cz,xj_cz,wc_mean] = support_set(wc,K);
sj_cz=sj_cz';xj_cz=xj_cz';wc_sj=wc_sj';wc_xj=wc_xj';
ee1=wasserstein(wc_mean,0.95,K,wc);
wc=wc';

[Mchp_fsbh,Mgb_fsbh,Mgrid_g_fsbh,Mgrid_s_fsbh,Mch_fsbh, Mdis_fsbh,Png_chp1_2_xj,Hng_chp1_2_xj,Hng_bo1_2_xj,Pgb1_2_xj,Pgs1_2_xj,Esp1_2_xj,Pd1_2_xj,Pc1_2_xj,...
    Png_chp1_2_sj,Hng_chp1_2_sj,Hng_bo1_2_sj,Pgb1_2_sj,Pgs1_2_sj,Esp1_2_sj,Pd1_2_sj,Pc1_2_sj,z1_dro,Pgb1_dro,Pgs1_dro,Pc1_dro,Pd1_dro,z_rqdro,mu,beta,...
    Png_chp1_dro,Hng_chp1_dro,Hng_bo1_dro,PL1d_dro,HL1d_dro,Ng1_dro,E1_dro,PL1_down_dro,PL1_cut_dro,HL1_up_dro,HL1_down_dro,y_sj,y_xj,y,GP,GH,HG,PG,Ng_chp1_dro,Ng_bo1_dro]...
    = SP_DRLQ(PV1,WT1,sj_cz,xj_cz,wc_xj,wc_sj,ee1,prgb,prgs,PL10,HL10,wc,GL10,GP_max,GH_max);

figure(1)
bar_data = [
    value(Pgb1_dro), value(Png_chp1_dro), value(Pd1_dro), value(WT1), value(PV1), value(-1*Pgs1_dro), value(-1*Pc1_dro), value(-1*PG)
];
b = bar(bar_data, 'stacked', 'BarWidth', 0.9);
b(1).FaceColor = [155/255, 190/255, 230/255];    
b(2).FaceColor = [175/255, 87/255, 99/255];    
b(3).FaceColor = [190/255, 180/255, 190/255]; 
b(4).FaceColor = [235/255, 210/255, 190/255];    
b(5).FaceColor = [255/255, 225/255, 190/255];    
b(6).FaceColor = [195/255, 120/255, 100/255];   
b(7).FaceColor = [150/255, 230/255, 220/255];   
b(8).FaceColor = [70/255, 190/255, 225/255];  
hold on;
plot(PL10, 'k', 'LineWidth', 1); 
xticks(4:4:192); 
xticklabels(mod(0:47, 24) + 1);
legend('电网购电', 'CHP出力', '储能放电', '海上风电出力', '光伏出力', '向电网售电', '储能充电', '电制冷机用电', '电负荷');
xlabel('时间/h');
ylabel('电功率/MW');
title('192小时电功率变化 (每15分钟调度)');
hold off;
figure(2)
b1 = bar([value(Hng_chp1_dro) value(Hng_bo1_dro) value(-1*HG)], 'stacked', 'BarWidth', 0.7);
hold on
b1(1).FaceColor = [0.8, 0.2, 0.3]; 
b1(2).FaceColor = [0.1, 0.6, 0.2]; 
b1(3).FaceColor = [0.9, 0.6, 0.1];

plot(HL10, 'r', 'LineWidth', 1); 

xticks(4:4:192);  
xticklabels(mod(0:47, 24) + 1);

legend('CHP供热', 'GB供热', 'AR吸热量', '热负荷');
xlabel('时间/h');
ylabel('电功率/W');
title('热网数据 (0-192 h)');
hold off
figure(3)
b3 = bar([value(Ng_bo1_dro) value(Ng_chp1_dro)], 'stacked', 'BarWidth', 0.7);
hold on
b3(1).FaceColor = [0.2, 0.6, 0.2]; 
b3(2).FaceColor = [0.1, 0.4, 0.8]; 
hold on
plot(Ng1_dro, 'g', 'LineWidth', 1); 
xticks(4:4:192);
xticklabels(mod(0:47, 24) + 1);
legend( 'GB消耗天然气', 'CHP消耗天然气','气网购气');
xlabel('时间/h');
ylabel('气功率/立方米');
title('气网数据 (0-192 h)');
hold off
figure(4)
b5 = bar([value(GH) value(GP) ], 'stacked', 'BarWidth', 0.7);
hold on
b5(1).FaceColor = [157/255, 215/255,157/255]; 
b5(2).FaceColor = [253/255,200/255, 151/255]; 
hold on
plot(GL10,'b','LineWidth', 1); 
xticks(4:4:192);  
xticklabels(mod(0:47, 24) + 1);

legend('AR制冷', 'EC制冷', '冷负荷');
xlabel('时间/h');
ylabel('冷功率/MW');
title('冷网数据 (0-192 h)');
hold off

figure(5)
plot(PL1d_dro,'g','linewidth',1.6)
hold on
plot(HL1d_dro,'r','linewidth',1.6)
hold on
plot(GL10,'b','linewidth',1.6)
hold on

xticks(4:4:192);  
xticklabels(mod(0:47, 24) + 1);
legend('电负荷','热负荷','冷负荷');
xlabel('时间/h');
ylabel('功率/kW');
figure(6)
plot(prgb,'g','linewidth',1.6)
hold on
xticks(4:4:192);  
xticklabels(mod(0:47, 24) + 1);
hold on
prgb22=readmatrix('仿真用_电价曲线模板.xlsx');
figure(7)
plot(prgb22,'g','linewidth',1.6)
hold on
xticks(4:4:192);  
xticklabels(mod(0:47, 24) + 1);
hold on
