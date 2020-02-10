# PyEMD

经验模态分解 (Empirical Mode Decomposition)
**Author: Aeo** 转载请注明出处

#### 算法背景
**经验模态分解**（Empirical Mode Decomposition，缩写EMD）是由[黄锷](https://baike.baidu.com/item/%E9%BB%84%E9%94%B7)（N. E. Huang）在美国国家宇航局与其他人于1998年创造性地提出的一种新型自适应信号时频处理方法，特别适用于非线性非平稳信号的分析处理。

---
#### 算法过程分析
- **筛选（Sifting）**
  1. **求极值点** 通过Find Peaks算法获取信号序列的全部`极大值`和`极小值`
  2. **拟合包络曲线** 通过信号序列的`极大值`和`极小值`组，经过`三次样条插值法`获得两条光滑的波峰/波谷拟合曲线，即信号的`上包络线`与`下包络线`
  3. **均值包络线** 将两条极值曲线平均获得`平均包络线`
  4. **中间信号** 原始信号减均值包络线，得到`中间信号`
  5. **判断本征模函数（IMF）** IMF需要符合两个条件：
        1）在整个数据段内，极值点的个数和过零点的个数必须相等或相差最多不能超过一个。
        2）在任意时刻，由局部极大值点形成的上包络线和由局部极小值点形成的下包络线的平均值为零，即上、下包络线相对于时间轴局部对称。
- **IMF 1** 获得的第一个满足IMF条件的中间信号即为原始信号的第一个本征模函数分量`IMF 1`
   （由原数据减去包络平均后的新数据，若还存在负的局部极大值和正的局部极小值，说明这还不是一个本征模函数，需要继续进行“筛选”。）
- 使用上述方法得到第一个IMF后，用原始信号减IMF1，作为新的原始信号，再通过上述的筛选分析，可以得到IMF2，以此类推，完成EMD分解。

---
#### 1. 求极值点
```
from scipy.signal import argrelextrema

# 通过Scipy的argrelextrema函数获取信号序列的极值点
data = np.random.random(100)
max_peaks = argrelextrema(data, np.greater)
min_peaks = argrelextrema(data, np.less)

# 绘制极值点图像
plt.figure(figsize = (18,6))
plt.plot(data)
plt.scatter(max_peaks, data[max_peaks], c='r', label='Max Peaks')
plt.scatter(min_peaks, data[min_peaks], c='b', label='Max Peaks')
plt.legend()
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
plt.title("Find Peaks")
```
![极大值与极小值点](https://upload-images.jianshu.io/upload_images/20053660-2b60bd9ab1feae2c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

---
#### 2. 拟合包络函数
这一步是EMD的核心步骤，也是分解出本征模函数IMFs的前提。
```
from scipy.signal import argrelextrema

data = np.random.random(100)-0.5
index = list(range(len(data)))

# 获取极值点
max_peaks = list(argrelextrema(data, np.greater)[0])
min_peaks = list(argrelextrema(data, np.less)[0])

# 将极值点拟合为曲线
ipo3_max = spi.splrep(max_peaks, data[max_peaks],k=3) #样本点导入，生成参数
iy3_max = spi.splev(index, ipo3_max) #根据观测点和样条参数，生成插值

ipo3_min = spi.splrep(min_peaks, data[min_peaks],k=3) #样本点导入，生成参数
iy3_min = spi.splev(index, ipo3_min) #根据观测点和样条参数，生成插值

# 计算平均包络线
iy3_mean = (iy3_max+iy3_min)/2

# 绘制图像
plt.figure(figsize = (18,6))
plt.plot(data, label='Original')
plt.plot(iy3_max, label='Maximun Peaks')
plt.plot(iy3_min, label='Minimun Peaks')
plt.plot(iy3_mean, label='Mean')
plt.legend()
plt.xlabel('time (s)')
plt.ylabel('microvolts (uV)')
plt.title("Cubic Spline Interpolation")
```
![三次样条插值法拟合包络函数](https://upload-images.jianshu.io/upload_images/20053660-6577b0463b7521fa.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

用原信号减去平均包络线即为所获得的新信号，若新信号中还存在负的局部极大值和正的局部极小值，说明这还不是一个本征模函数，需要继续进行“筛选”。
![第一次筛选后的新信号与原始信号对比](https://upload-images.jianshu.io/upload_images/20053660-60f7e4843afa8cb2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


---
#### 获取本征模函数（IMF）
```
def sifting(data):
    index = list(range(len(data)))

    max_peaks = list(argrelextrema(data, np.greater)[0])
    min_peaks = list(argrelextrema(data, np.less)[0])

    ipo3_max = spi.splrep(max_peaks, data[max_peaks],k=3) #样本点导入，生成参数
    iy3_max = spi.splev(index, ipo3_max) #根据观测点和样条参数，生成插值

    ipo3_min = spi.splrep(min_peaks, data[min_peaks],k=3) #样本点导入，生成参数
    iy3_min = spi.splev(index, ipo3_min) #根据观测点和样条参数，生成插值

    iy3_mean = (iy3_max+iy3_min)/2
    return data-iy3_mean


def hasPeaks(data):
    max_peaks = list(argrelextrema(data, np.greater)[0])
    min_peaks = list(argrelextrema(data, np.less)[0])
    
    if len(max_peaks)>3 and len(min_peaks)>3:
        return True
    else:
        return False


# 判断IMFs
def isIMFs(data):
    max_peaks = list(argrelextrema(data, np.greater)[0])
    min_peaks = list(argrelextrema(data, np.less)[0])
    
    if min(data[max_peaks]) < 0 or max(data[min_peaks])>0:
        return False
    else:
        return True

    
def getIMFs(data):
    while(not isIMFs(data)):
        data = sifting(data)
    return data


# EMD函数
def EMD(data):
    IMFs = []
    while hasPeaks(data):
        data_imf = getIMFs(data)
        data = data-data_imf
        IMFs.append(data_imf)
    return IMFs


# 绘制对比图
data = np.random.random(1000)-0.5
IMFs = EMD(data)
n = len(IMFs)+1

# 原始信号
plt.figure(figsize = (18,15))
plt.subplot(n, 1, 1)
plt.plot(data, label='Origin')
plt.title("Origin ")

# 若干条IMFs曲线
for i in range(0,len(IMFs)):
    plt.subplot(n, 1, i+2)
    plt.plot(IMFs[i])
    plt.ylabel('Amplitude')
    plt.title("IMFs "+str(i+1))

plt.legend()
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
```
![迭代分解结果](https://upload-images.jianshu.io/upload_images/20053660-6f59ed075803a675.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

至此就完成了**经验模态分解EMD**的基本过程。

---
#### Reference
[【百度百科】经验模态分解](https://baike.baidu.com/item/%E7%BB%8F%E9%AA%8C%E6%A8%A1%E6%80%81%E5%88%86%E8%A7%A3/1238591?fr=aladdin)

[【知乎】这篇文章能让你明白经验模态分解（EMD）——基础理论篇](https://zhuanlan.zhihu.com/p/40005057)

[【知乎】这篇文章能让你明白经验模态分解（EMD）——IMF的物理含义](https://zhuanlan.zhihu.com/p/44833026)


[【CSDN】Python实现线性插值和三次样条插值](https://blog.csdn.net/weixin_41799019/article/details/97629116)
