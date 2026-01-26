function[Mchp_fsbh,Mgb_fsbh,Mgrid_g_fsbh,Mgrid_s_fsbh,Mch_fsbh, Mdis_fsbh,Png_chp1_2_xj,Hng_chp1_2_xj,Hng_bo1_2_xj,Pgb1_2_xj,Pgs1_2_xj,Esp1_2_xj,Pd1_2_xj,Pc1_2_xj,...
    Png_chp1_2_sj,Hng_chp1_2_sj,Hng_bo1_2_sj,Pgb1_2_sj,Pgs1_2_sj,Esp1_2_sj,Pd1_2_sj,Pc1_2_sj,z1_dro,Pgb1_dro,Pgs1_dro,Pc1_dro,Pd1_dro,z_rqdro,mu,beta,...
    Png_chp1_dro,Hng_chp1_dro,Hng_bo1_dro,PL1d_dro,HL1d_dro,Ng1_dro,E1_dro,PL1_down_dro,PL1_cut_dro,HL1_up_dro,HL1_down_dro,y_sj,y_xj,y,GP,GH,HG,PG,Ng_chp1_dro,Ng_bo1_dro]...
    = SP_DRLQ(PV1,WT1,sj_cz,xj_cz,wc_xj,wc_sj,ee1,prgb,prgs,PL10,HL10,wc,GL10,GP_max,GH_max)
Png_chp1_max=1000;
Hng_bo1_max=700;
Esp1_max=50;
pd1_shuzhi=30;
pc1_shuzhi=30;
eff_ch=0.95;
eff_dch=0.95;
D1=0.8;
D2=0.85;
LHV=0.00978;
chp_up=0.1;
gb_up=0.05;
dr1=0;
T=192;
COP_EC=3.5; 
COP_AC=1.2;

Pgb1= sdpvar(T,1);
Pgs1= sdpvar(T,1);
Ng1= sdpvar(T,1);
ng_price=4;
PL1d=sdpvar(T,1);
HL1d=sdpvar(T,1);

Esp1 = sdpvar(T,1);
Pc1 = sdpvar(T,1);
Pd1 = sdpvar(T,1);
binC1 = binvar(T,1);

Ng_chp1=sdpvar(T,1);
Png_chp1= sdpvar(T,1);
yita_chpd=0.4;
yita_chph=0.4;
yita_chuanre=0.6;
Hng_chp1 = sdpvar(T,1);
Binchp1_ng=binvar(T,1); 
Hng_bo1 = sdpvar(T,1);
price_gb=0.01;
PL1_cut=sdpvar(T,1);
PL1_up=sdpvar(T,1);
PL1_down=sdpvar(T,1);
PL1_cv=sdpvar(T,1);
HL1_cut=sdpvar(T,1);
HL1_up=sdpvar(T,1);
HL1_down=sdpvar(T,1);
HL1_cv=sdpvar(T,1);
B1=binvar(T,1);
B2=binvar(T,1);
E1 = sdpvar(T,1);
a1=binvar(T,1);
a2=binvar(T,1);
a3=binvar(T,1);
a4=binvar(T,1);
a5=binvar(T,1);
a6=binvar(T,1);
a7=binvar(T,1);
a8=binvar(T,1);
a9=binvar(T,1);
M=10000;

GP=sdpvar(T,1);
GH=sdpvar(T,1);
PG=sdpvar(T,1);
HG=sdpvar(T,1);
 mu=sdpvar(1);
 beta= sdpvar(K,1,'full');
Mchp_fsbh=sdpvar(T,1); 
Mgb_fsbh=sdpvar(T,1);  
Mgrid_g_fsbh=sdpvar(T,1);  
Mgrid_s_fsbh=sdpvar(T,1);  
Mch_fsbh=sdpvar(T,1); 
Mdis_fsbh=sdpvar(T,1);  
y=sdpvar(K,1);
Png_chp1_2_sj=sdpvar(T,1); 
Hng_chp1_2_sj=sdpvar(T,1);
Hng_bo1_2_sj=sdpvar(T,1);
Pgb1_2_sj=sdpvar(T,1);
Esp1_2_sj=sdpvar(T,1);
Pd1_2_sj=sdpvar(T,1);
Pc1_2_sj=sdpvar(T,1);
b2_sj=binvar(T,1);
b3_sj=binvar(T,1);
b4_sj=binvar(T,1);
Png_chp1_2_xj=sdpvar(T,1); 
Hng_bo1_2_xj=sdpvar(T,1);
Pgb1_2_xj=sdpvar(T,1);
Pgs1_2_xj=sdpvar(T,1);
Esp1_2_xj=sdpvar(T,1);
Pd1_2_xj=sdpvar(T,1);
Pc1_2_xj=sdpvar(T,1);
b4_xj=sdpvar(T,1);
b2_xj=binvar(T,1);
b3_xj=binvar(T,1);
cons = [];
for i = 1:K
    ycwc = wc(:,i);

    Ng1_2 = -Mchp_fsbh .* ycwc / (LHV*yita_chpd) + (-Mgb_fsbh .* ycwc) / (LHV*0.85);

    y(i) = sum( ...
        prgb .* (-Mgrid_g_fsbh .* ycwc) * 250 ...
        - prgs .* (Mgrid_s_fsbh .* ycwc) * 250 ...
        + Ng1_2 * ng_price*0.25 ...
        + 0.07 * (Mch_fsbh .* ycwc - Mdis_fsbh .* ycwc) * 250 ...
    );
end
Ng1_2_sj = -Mchp_fsbh .* wc_sj / (LHV*yita_chpd) + (-Mgb_fsbh .* wc_sj) / (LHV*0.85);

y_sj = sum( ...
    prgb .* (-Mgrid_g_fsbh .* wc_sj) * 250 ...
    - prgs .* (Mgrid_s_fsbh .* wc_sj) * 250 ...
    + Ng1_2_sj * ng_price*0.25 ...
    + 0.07 * (Mch_fsbh .* wc_sj - Mdis_fsbh .* wc_sj) * 250 ...
);
Ng1_2_xj = -Mchp_fsbh .* wc_xj / (LHV*yita_chpd) + (-Mgb_fsbh .* wc_xj) / (LHV*0.85);

y_xj = sum( ...
    prgb .* (-Mgrid_g_fsbh .* wc_xj) * 250 ...
    - prgs .* (Mgrid_s_fsbh .* wc_xj) * 250 ...
    + Ng1_2_xj * ng_price*0.25 ...
    + 0.07 * (Mch_fsbh .* wc_xj - Mdis_fsbh .* wc_xj) * 250 ...
);
for j=1:K
cons = [cons,beta(j)>=y(j)];
cons = [cons,beta(j)>=y_sj-mu.*ones(1,192)*sj_cz(:,j)];
cons = [cons,beta(j)>=y_xj-mu.*ones(1,192)*xj_cz(:,j)];
end
cons = [cons,mu>=0];
cons=[cons,0<=Png_chp1];
cons=[cons,-chp_up*Png_chp1_max<=Png_chp1(2:T)-Png_chp1(1:T-1)<=chp_up*Png_chp1_max];
cons=[cons,Ng_chp1*LHV*yita_chpd==Png_chp1];
cons=[cons,Ng_chp1*LHV*yita_chph*yita_chuanre==Hng_chp1];
cons=[cons,Ng_bo1*LHV*0.85==Hng_bo1];
cons=[cons,-gb_up*Hng_bo1_max<=Hng_bo1(2:T)-Hng_bo1(1:T-1)<=gb_up*Hng_bo1_max];
cons=[cons,Hng_bo1==0];
cons=[cons,Hng_bo1<=Hng_bo1_max];
cons=[cons, Esp1>=(1-D1)*Esp1_max*ones(T,1)];
cons=[cons, Esp1<=0.95*Esp1_max*ones(T,1)];
cons = [cons, Esp1(1) == (1-D1)*Esp1_max + Pc1(1)*eff_ch*0.25];
cons=[cons, Pd1(1)==0];
cons = [cons, Esp1(2:192) == Esp1(1:191) + (Pc1(2:192)*eff_ch - Pd1(2:192)/eff_dch)*0.25];
cons = [cons, Esp1(192)==(1-D1)*Esp1_max*ones(T,1)];
cons=[cons,Pd1>=0];
cons=[cons,Pc1>=0];
cons = [cons, Pd1<=pd1_shuzhi];
cons = [cons, Pc1<=pc1_shuzhi];
cons=[cons,Pc1<=binC1*M];
cons=[cons,Pc1(141:153)==0];
cons=[cons,PL1_cut<=dr1*PL10];
cons=[cons,PL1_cut>=0];
cons=[cons,sum(PL1_up)-sum(PL1_down)==0];
cons=[cons,PL1_up<=dr1*PL10];
cons=[cons,PL1_up>=0];
cons=[cons,PL1_down<=dr1*PL10];
cons=[cons,PL1_down>=0];
cons=[cons,PL1_up<=(1-B1)*M];
cons=[cons,PL1_down<=B1*M];
cons=[cons,PL1_cv<=dr1*PL10];
cons=[cons,PL1_cv>=-dr1*PL10];
cons=[cons,HL1_cv==2*PL1_cv];
cons=[cons,HL1_cut<=dr1*HL10];%热切负荷
cons=[cons,HL1_cut>=0];  
cons=[cons,sum(HL1_up)-sum(HL1_down)==0];
cons=[cons,HL1_up<=dr1*HL10];
cons=[cons,HL1_up>=0];
cons=[cons,HL1_down<=dr1*HL10]; 
cons=[cons,HL1_down>=0];
cons=[cons,HL1_up<=(1-B2)*M];
cons=[cons,HL1_down<=B2*M];
cons=[cons,PL1d==PL10-PL1_cut+PL1_up-PL1_down+PL1_cv];%
cons=[cons,HL1d==HL10-HL1_cut+HL1_up-HL1_down-HL1_cv];%
cons=[cons,PL1d>=0];
cons=[cons,HL1d>=0];
cons=[cons,GP==PG*COP_EC]; 
cons=[cons,GP>=0];
cons=[cons,GP<=GP_max];
cons=[cons,GH==HG*COP_AC]; 
cons=[cons,GH>=0];
cons=[cons,GH<=GH_max];
cons=[cons,Pgb1+WT1+PV1+Png_chp1+Pd1==Pgs1+PL10+Pc1+PG];
cons=[cons,Hng_chp1+Hng_bo1-HL10-HG==0];
cons=[cons,Ng1==Ng_bo1+Ng_chp1];
cons=[cons,GP+GH-GL10==0];
cons=[cons,Pgb1>=0]
cons=[cons,Pgs1==0];
cons=[cons,Pgb1<=(1-a1)*M];
cons=[cons,Pgs1<=a1*M];

cons=[cons,Mchp_fsbh>=-1];cons=[cons,Mchp_fsbh<=1];cons=[cons,Mgrid_g_fsbh>=-1];cons=[cons,Mgrid_g_fsbh<=1];
cons=[cons,Mgrid_s_fsbh>=-1];cons=[cons,Mgrid_s_fsbh<=1]; %-1<上网电量调整系数<1
cons=[cons,Mch_fsbh>=-1];cons=[cons,Mch_fsbh<=1];cons=[cons,Mdis_fsbh>=-1];cons=[cons,Mdis_fsbh<=1];cons=[cons,Mgb_fsbh>=-1];cons=[cons,Mgb_fsbh<=1];

cons=[cons,Png_chp1_2_sj==Png_chp1-Mchp_fsbh.*wc_sj];
cons=[cons,-chp_up*Png_chp1_max<=Png_chp1_2_sj(2:T)-Png_chp1_2_sj(1:T-1)<=chp_up*Png_chp1_max];
cons=[cons,0<=Png_chp1_2_sj];
cons=[cons,Hng_chp1_2_sj==Png_chp1_2_sj*3*yita_chuanre*yita_chph];%CHP机组放热功率上界=CHP机组发电功率上界*3*0.47

cons=[cons,Hng_bo1_2_sj==Hng_bo1-Mgb_fsbh.*wc_sj];%燃气锅炉放热功率上界=GB-调整系数*误差
cons=[cons,-gb_up*Hng_bo1_max<=Hng_bo1_2_sj(2:T)-Hng_bo1_2_sj(1:T-1)<=gb_up*Hng_bo1_max];
cons=[cons,Hng_bo1_2_sj>=0];
cons=[cons,Hng_bo1_2_sj<=Hng_bo1_max];

cons=[cons,Pgb1_2_sj==Pgb1-Mgrid_g_fsbh.*wc_sj];
cons=[cons,Pgb1_2_sj>=0];
cons=[cons,Pgs1_2_sj==Pgs1+Mgrid_s_fsbh.*wc_sj];
cons=[cons,Pgs1_2_sj>=0];
cons=[cons,Pgs1_2_sj<=(1-b3_sj)*M];
cons=[cons,Pgb1_2_sj<=(b3_sj)*M];

cons=[cons, Pd1_2_sj(1)==0];
cons=[cons, Mdis_fsbh(1)==0];
cons=[cons, Pd1_2_sj(2:T)==Pd1(2:T)-Mdis_fsbh(2:T).*wc_sj(2:T)];
cons=[cons, Pc1_2_sj==Pc1+Mch_fsbh.*wc_sj];
cons=[cons, Esp1_2_sj>=(1-D1)*Esp1_max*ones(T,1)];
cons=[cons, Esp1_2_sj<=0.95*Esp1_max*ones(T,1)];
cons=[cons, Esp1_2_sj(1)== (1-D1)*Esp1_max + Pc1_2_sj(1)*eff_ch*0.25];
cons=[cons, Esp1_2_sj(2:192)==Esp1_2_sj(1:191)+(Pc1_2_sj(2:192)*eff_ch - Pd1_2_sj(2:192)/eff_dch)*0.25];
cons=[cons, Esp1_2_sj(192)==(1-D1)*Esp1_max];
cons=[cons, Pd1_2_sj>=0];
cons=[cons, Pc1_2_sj>=0];
cons=[cons, Pd1_2_sj<=pd1_shuzhi];
cons=[cons, Pc1_2_sj<=pc1_shuzhi];
cons=[cons, Pd1_2_sj<=(1-b2_sj)*M];
cons=[cons, Pc1_2_sj<=b2_sj*M];
cons=[cons,Pc1_2_sj(141:153)==0];

cons=[cons,-wc_sj-Mchp_fsbh.*wc_sj-Mgrid_g_fsbh.*wc_sj-Mgrid_s_fsbh.*wc_sj-Mdis_fsbh.*wc_sj-Mch_fsbh.*wc_sj==0];
cons=[cons,-Mchp_fsbh.*wc_sj*yita_chuanre*yita_chph*3-Mgb_fsbh.*wc_sj==0];%这里有购热系数

cons=[cons,Png_chp1_2_xj==Png_chp1-Mchp_fsbh.*wc_xj];%
cons=[cons,-chp_up*Png_chp1_max<=Png_chp1_2_xj(2:T)-Png_chp1_2_xj(1:T-1)<=chp_up*Png_chp1_max];
cons=[cons,0<=Png_chp1_2_xj];
cons=[cons,Hng_chp1_2_xj==Png_chp1_2_xj*3*yita_chuanre*yita_chph];

cons=[cons,Hng_bo1_2_xj==Hng_bo1-Mgb_fsbh.*wc_xj];
cons=[cons,-gb_up*Hng_bo1_max<=Hng_bo1_2_xj(2:T)-Hng_bo1_2_xj(1:T-1)<=gb_up*Hng_bo1_max];
cons=[cons,Hng_bo1_2_xj>=0];
cons=[cons,Hng_bo1_2_xj<=Hng_bo1_max];

cons=[cons,Pgb1_2_xj==Pgb1-Mgrid_g_fsbh.*wc_xj];
cons=[cons,Pgb1_2_xj>=0];
cons=[cons,Pgs1_2_xj==Pgs1+Mgrid_s_fsbh.*wc_xj];
cons=[cons,Pgs1_2_xj>=0];
cons=[cons,Pgs1_2_xj<=(1-b3_xj)*M];
cons=[cons,Pgb1_2_xj<=(b3_xj)*M];
cons=[cons, Pd1_2_xj(1)==0];
cons=[cons, Mdis_fsbh(1)==0];
cons=[cons, Pd1_2_xj(2:T)==Pd1(2:T)-Mdis_fsbh(2:T).*wc_xj(2:T)];
cons=[cons, Pc1_2_xj==Pc1+Mch_fsbh.*wc_xj];
cons=[cons, Esp1_2_xj>=(1-D1)*Esp1_max*ones(T,1)];
cons=[cons, Esp1_2_xj<=0.95*Esp1_max*ones(T,1)];
cons=[cons, Esp1_2_xj(1)== (1-D1)*Esp1_max + Pc1_2_xj(1)*eff_ch*0.25];
cons=[cons, Esp1_2_xj(2:192)==Esp1_2_xj(1:191)+(Pc1_2_xj(2:192)*eff_ch - Pd1_2_xj(2:192)/eff_dch)*0.25];
cons=[cons, Esp1_2_xj(192)==(1-D1)*Esp1_max];
cons=[cons, Pd1_2_xj>=0];
cons=[cons, Pc1_2_xj>=0];
cons=[cons, Pd1_2_xj<=pd1_shuzhi];
cons=[cons, Pc1_2_xj<=pc1_shuzhi];
cons=[cons, Pd1_2_xj<=(1-b2_xj)*M];
cons=[cons, Pc1_2_xj<=b2_xj*M];
cons=[cons,Pc1_2_xj(141:153)==0];

cons=[cons,-wc_xj-Mchp_fsbh.*wc_xj-Mgrid_g_fsbh.*wc_xj-Mgrid_s_fsbh.*wc_xj-Mdis_fsbh.*wc_xj-Mch_fsbh.*wc_xj==0];
cons=[cons,-Mchp_fsbh.*wc_xj*yita_chuanre*yita_chph*3-Mgb_fsbh.*wc_xj==0];

rq = sum( ...
    (prgb .* Pgb1 * 250 - prgs .* Pgs1 * 250) ... 
    + (Ng1 * ng_price) ... 
    + (price_chp * Png_chp1 * 250) ... 
    + (price_gb * Hng_bo1 * 250) ... 
    + (0.05 * (Pd1 + Pc1) * 250) ... 
    + (0.07 * (PL1_up + PL1_down) * 250) ... 
    + (0.14 * PL1_cut * 250) ... 
    + (0.035 * (HL1_up + HL1_down) * 250) ... 
    + (0.07 * HL1_cut * 250) ... 
    + (0.05 * GP * 250) ... 
    + (0.04 * GH * 250) ... 
);
obj=rq+(sum(beta))/K+mu*ee1;

fprintf('solver')
options=sdpsettings('solver','gurobi'); 
options.debug=1;
options.showprogress=1;
optimize(cons,obj,options);
beta=value(beta);mu=value(mu);
z1_dro=value(obj);
Mchp_fsbh=value(Mchp_fsbh);Mgb_fsbh=value(Mgb_fsbh);Mgrid_g_fsbh=value(Mgrid_g_fsbh);Mgrid_s_fsbh=double(Mgrid_s_fsbh); Mch_fsbh=double(Mch_fsbh); Mdis_fsbh=double(Mdis_fsbh);
Pgb1_2_xj=value(Pgb1_2_xj);Pgb1_2_sj=value(Pgb1_2_sj);Png_chp1_2_xj=value(Png_chp1_2_xj);Png_chp1_2_sj=value(Png_chp1_2_sj);Hng_bo1_2_xj=value(Hng_bo1_2_xj);Hng_bo1_2_sj=value(Hng_bo1_2_sj);Pgs1_2_sj=value(Pgs1_2_sj);Pgs1_2_xj=value(Pgs1_2_xj);
Hng_chp1_2_sj=value(Hng_chp1_2_sj);Hng_chp1_2_xj=value(Hng_chp1_2_xj);
Esp1_2_sj=value(Esp1_2_sj);Esp1_2_xj=value(Esp1_2_xj);Pd1_2_sj=value(Pd1_2_sj);Pd1_2_xj=value(Pd1_2_xj);Pc1_2_sj=value(Pc1_2_sj);Pc1_2_xj=value(Pc1_2_xj);
Pgb1_dro=value(Pgb1);Pgs1_dro=value(Pgs1);Png_chp1_dro=value(Png_chp1);Hng_chp1_dro=value(Hng_chp1);Hng_bo1_dro=value(Hng_bo1);
PL1d_dro=double(PL1d); HL1d_dro=double(HL1d);Ng1_dro=value(Ng1);E1_dro=value(E1);z_rqdro=value(rq);
Pc1_dro=double(Pc1); Pd1_dro=double(Pd1);PL1_down_dro=double(PL1_down);PL1_cut_dro=double(PL1_cut); 
HL1_up_dro=double(HL1_up);HL1_down_dro=double(HL1_down);HL1_cut_dro=double(HL1_cut);PL1_up_dro=double(PL1_up);
PL1_cv_dro=double(PL1_cv);HL1_cv_dro=double(HL1_cv);y=value(y);y_sj=value(y_sj);y_xj=value(y_xj);
GP=value(GP);GH=value(GH);PG=value(PG);HG=value(HG);Ng_chp1_dro=value(Ng_chp1);Ng_bo1_dro=value(Ng_bo1);
