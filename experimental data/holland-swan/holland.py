from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import datetime as dt
import netCDF4 as nc4
plt.rcParams['font.sans-serif']=['Times New Roman']#设置绘图字体为 Times New Roman，确保生成的图像中文或英文的字体统一。
#读取nc文件
file2 = Dataset('.\\bajiao_era.nc')   #修改1，era文件名
lon = file2['longitude'][:]
lat = file2['latitude'][:]

u10 = file2['u10'][:]
v10 = file2['v10'][:]
p10 = file2['sp'][:]#地表气压
time0 =file2['valid_time'][:]#时间信息
time2=[]
tstart = dt.datetime(1970,1,1,0)
for i in time0:
        time2.append(tstart+dt.timedelta(hours=int(i/3600)))#将时间变量从相对于 1970 年的秒数转化为人类可读的日期时间格式。
print(time2)
Lon,Lat = np.meshgrid(lon,lat)        #
data = np.loadtxt('.\\2022芭蕉.txt',delimiter='\t',usecols=(0,1,2,3))   #修改2，txt文件名第 0 列：台风中心纬度。第 1 列：台风中心经度。第 2 列：台风中心气压。第 3 列：最大风速。
print(data)
T1 = range(0,len(data[:,0])*6,6);  #
T2 = range(0,(len(data[:,0])-1)*6+1,1);

lat2 = np.interp(T2, T1,data[:,0]);   #插值后的纬度
long2 =  np.interp(T2,T1, data[:,1]);#插值后的经度
p =  np.interp(T2,T1, data[:,2]);#插值后的气压
velc =  np.interp(T2,T1, data[:,3]);#插值后的风速

vel = np.zeros((len(lat),len(lon)))#初始化各变量为全零矩阵，用于存储计算结果
u = np.zeros((len(lat),len(lon)))
v = np.zeros((len(lat),len(lon)))
p0 = np.zeros((len(lat),len(lon)))
uout = np.zeros((len(lat2),len(lat),len(lon)))
vout = np.zeros((len(lat2),len(lat),len(lon)))
pout = np.zeros((len(lat2),len(lat),len(lon)))
vfmax= []
velmaxcal = []
for i in range(len(T2)-1):
        fai = lat2[i] / 180 * np.pi
        vfmax.append(np.sqrt((lat2[i]-lat2[i+1])**2+((long2[i]-long2[i+1])*np.cos(fai))**2)*111320/3600/6)
vfmax.append(vfmax[-1])#计算每个时间步台风路径点的最大风速（vfmax），考虑了台风路径的移动速度。

num=0
for i in range(0,len(T2),1):
        print(i)
        Pc = p[num]*100
        fai = lat2[i]/180*np.pi
        Vfmax = vfmax[num]
        rou = 1.29
        omiga = 7.292e-5
        f = 2 * omiga * np.sin(fai)
        Rmax=28.52*np.tanh(0.0873*(lat2[i]-28))+12.22*np.exp((p[i]-1013.2)/33.86)+0.2*Vfmax/1000+37.2
        #Rmax = 51.6 * np.exp(-0.022 * Vfmax / 1000 + 0.028 * fai)
        #Rmax = np.exp(2.636-0.00005086*(p[i]-1013.2)**2+0.0394899*fai/np.pi*180)
        #Rmax = np.exp(3.94-0.0223*velc[i]+0.0281*fai/np.pi*180)
        #Rmax = np.exp(2.0633+0.0182*(p[i]-1013.2)-0.00019008*(p[i]-1013.2)**2+0.0007336*(fai/np.pi*180)**2)
        print('ramx',Rmax)
        Rmax *= 1000
        B = 1.5+(980-Pc/100)/120#
        #B = (velc[i]**2+Rmax/1000*f*velc[i])*rou*2.713/(1010-p[i])/100
        #B = (velc[i] ** 2) * rou * 2.713 / (1010 - p[i]) / 100
        #B = 1.881-0.006*Rmax/1000-0.013*fai
        # B=0.25+0.3*np.log(1010-Pc/100)
        # B=1.38-0.00184*(1010-Pc/100)+0.00309*Rmax/1000
        vel[:,:] = 0
        for j in range(len(lat)):
                for k in range(len(lon)):
                        y = lat2[i]
                        x = long2[i]
                        r = np.sqrt((lat[j]-y)**2+((lon[k]-x)*np.cos(fai))**2)*111320
                        vel[j,k]=(B/1.2*(Rmax/r)**B*((1010-p[i])*100)*np.exp(-(Rmax/r)**B)+(r*f/2)**2)**0.5-r*f/2
                        p0[j, k] = (p[i] + (1013.2 - p[i]) * np.exp(-(Rmax / r) ** B))*100
                        if r<=1:
                                vel[j,k]=0
                        if lon[k]-long2[i] !=0:
                                theta = np.arctan(abs((lat[j]-lat2[i]))/abs((lon[k]-long2[i])))
                        else:
                                theta = np.pi*0.5
                        if (lat[j]-lat2[i])>=0 and (lon[k]-long2[i])>0:
                                u[j, k] = vel[j, k] * (-np.sin(theta))
                                v[j,k]= vel[j,k]*np.cos(theta)
                        elif (lat[j]-lat2[i])>=0 and (lon[k]-long2[i])<=0:
                                u[j, k] = vel[j, k] * (-np.sin(theta))
                                v[j, k] = -vel[j, k] * np.cos(theta)
                        elif (lat[j]-lat2[i])<0 and (lon[k]-long2[i])<=0:
                                u[j, k] = vel[j, k] * (np.sin(theta))
                                v[j, k] = -vel[j, k] * np.cos(theta)
                        elif (lat[j] - lat2[i]) < 0 and (lon[k] -long2[i]) >= 0:
                                u[j, k] = vel[j, k] * (np.sin(theta))
                                v[j, k] = vel[j, k] * np.cos(theta)
                        C=r/9/Rmax
                        if C<=5/9:
                                e= C**4/(1+C**4)
                                u[j, k] = u[j, k] * (1 - e) + u10[num, j, k] * e#
                                v[j, k] = v[j, k] * (1 - e) + v10[num, j, k] * e#
                                p0[j, k] = p0[j, k] * (1 - e) + p10[num, j, k] * e
                        else:
                                u[j, k]=u10[num,j,k]
                                v[j, k] = v10[num,j,k]
                                p0[j, k] = p10[num, j, k]

        vel=np.sqrt(u**2+v**2)
        uout[i]=u
        vout[i]=v
        pout[i]=p10[i]
        map = Basemap(llcrnrlon=lon[0], llcrnrlat=lat[-1], urcrnrlon=lon[-1], urcrnrlat=lat[0])  # 建立map对象
        lon = np.array(lon)
        lat = np.array(lat)  # 转为np
        if i%6==0:#可视化风场
                fig, ax = plt.subplots(figsize=(10, 7))  # 建立fig对象
                map.readshapefile(".\\china-shapefiles-master\\shapefiles\\ne_50m_admin_0_countries", 'china',
                                drawbounds=True, linewidth=1.5)  # 打开就又岸线
                fig = map.contourf(Lon, Lat, vel, np.linspace(10, 60, 51), cmap='plasma')
                # fig = map.contourf(Lon, Lat, vel, np.linspace(0,60,61),cmap='rainbow')  # fig填充map网格,cmap为色条种类
                font = {
                        'size': 20,
                }
                cbar = map.colorbar(fig, location='right', pad='10%')  # 图例
                cbar.ax.tick_params(labelsize=20)
                cbar.set_label('m/s', fontdict=font)
                plt.quiver(lon, lat, u, v, width=0.001, scale=150, scale_units='inches',
                       color='k', minshaft=2)
                plt.scatter(data[:, 1], data[:, 0], color='red', edgecolor='black', s=50)
                plt.plot(data[:, 1], data[:, 0], color='red', linewidth=2, linestyle='-', label='Typhoon Track')
                #plt.legend(fontsize=12)
                #plt.scatter(data[:, 1], data[:, 0],color='w')
                map.drawparallels(np.arange(lat[-1], lat[0] + 0.1, 2), labels=[1, 0, 0, 0], fontsize=20, alpha=1.0,
                                 color='none')  # 纬度线
                map.drawmeridians(np.arange(lon[0], lon[-1] + 0.1, 5), labels=[0, 0, 0, 1], fontsize=20, alpha=1.0,
                                  color='none')  # 经度线
                plt.title(str(time2[num]), fontsize=20, )  # 图名，英文
                print(np.max(vel))
                velmaxcal.append(np.max(vel))
                plt.savefig('vel_holland_era' + str(i) + '.jpg')
                ax.clear()
                #plt.show()

        num += 1
name = '2022bajiao_holland_era.nc'   #修改3，输出文件名
print(name)
ncfile = nc4.Dataset(name, 'w', format='NETCDF4')#保存NetCDF文件
ncfile.createDimension('lat', len(lat))
ncfile.createDimension('long', len(lon))
ncfile.createDimension('time', len(lat2))
time = ncfile.createVariable('time', np.float64, ('time',))
long = ncfile.createVariable('long', np.float64, ('long',))
LAT = ncfile.createVariable('lat', np.float64, ('lat',))
u10 = ncfile.createVariable('u10', np.float64, ('time', 'lat', 'long',))
v10 = ncfile.createVariable('v10', np.float64, ('time','lat', 'long', ))
sp = ncfile.createVariable('sp', np.float64, ('time','lat', 'long', ))
LAT.units = 'degree_north'
long.units = 'degrees_east'
LAT.long_name = 'latitude'
long.long_name ='longitude'
time.units = 'hours since 1900-01-01'
time.long_name = 'time'
long[:] = lon[:]
LAT[:] = lat[:]
time[:]=time0[:len(lat2)]
print(np.shape(u10))
print(np.shape(uout))
u10[...] = uout
v10[...] = vout
sp[...] = pout
ncfile.close()
print(velmaxcal)
np.savetxt('./maxvel_holland叠加风场.dat',np.array(velmaxcal))#保存最大风速
