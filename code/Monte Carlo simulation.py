import numpy as np
import scipy.stats as stats
from scipy.stats import beta
import time
start_time = time.time()  # 获取当前时间

# 步骤1：设备参数初始化
# 制造系统参数
R_max = 16.25  # 最大生产率
D = 10  # 固定需求率
α_X = 0.06  # 加速退化伽马过程的形状参数
β_X = 0.64  # 加速退化伽马过程的尺度参数
L = 7.5  # 进行CM的劣化阈值
X = 0  # 系统劣化水平
I = 0  # 初始库存
# 维护参数
u = 1  # β分布的参数
v = 5  # β分布的参数
g = 1.045  # 几何过程比率
'''T_CM = 0
T_PM = 0
T_OH = 0'''
T_Total = 0
# 质量相关参数
α = 1  # 浴盆曲线的尺度参数
β = 0.00025  # 浴盆曲线的尺度参数
θ = 0.05  # 浴盆曲线的形状参数
γ = 0.01  # 浴盆曲线的形状参数
λ = 1.75  # 浴盆曲线的加速度参数
p_0 = 0.02  # 缺陷率的初值
p_1 = 0.28  # 缺陷率的峰值
m_0 = 0.02  # 返修失败率的初值
m_1 = 0.28  # 返修失败率的峰值
AOQL = 0.0125  # 平均输出质量限制
AOQ = 0
N_def = 0 # 当前生产的缺陷产品数
N_rep = 0 # 当前生产的需要返工的产品数
N_ins = 0 # 到目前为止的检查次数
N_tcb = 0 # 到目前为止发送给顾客的坏品数
N_tc = 0 # 到目前为止发送给顾客的产品数
N_tcb_new = 0 # 到目前为止发送给顾客的坏品数
# 决策变量
Z = 320
T_E = 30
ξ = 5.4
Is = 16
# 成本参数
c_set = 800
c_PM = 600  # PM成本
c_CM = 2000  # CM成本
c_OH = 1800  # 大修成本
c_rep = 30  # 缺陷产品的维修成本
c_rej = 50  # 返工缺陷产品的拒绝成本
c_h = 1  # 单位时间产品的持有成本
c_s = 120  # 单位时间产品的短缺成本
c_ins = 12.5  # 单次检查成本
c_q = 150  # 出售给客户不良商品的单位成本
C_Total = 0
C_h = 0
# 仿真参数
t = 0
T = 10000  # 仿真总时间
n_1 = 0  # Scenario1循环出现的次数
n_2 = 0  # Scenario2循环出现的次数
n_3 = 0  # Scenario3循环出现的次数
n_4 = 0  # Scenario4循环出现的次数
n_5 = 0  # Scenario5循环出现的次数
n_6 = 0  # Scenario6循环出现的次数
n_oldsystem = 0
k = 0
set = 0  # 生产运行的设置成本
PM = 0  # PM成本
CM = 0  # CM成本
OH = 0  # 大修成本
rep = 0  # 缺陷产品的维修成本
rej = 0  # 返工缺陷产品的拒绝成本
h = 0  # 单位时间产品的持有成本
s = 0  # 单位时间产品的短缺成本
ins = 0  # 单次检查成本
q = 0
epsilon0 = 1e-5  # 定义一个非常小的正数，用来避免除以零
epsilon1 = 1e-5  # 定义一个非常小的正数，用来避免除以零
# 步骤2：模拟系统的退化过程
# 2.1
while T_Total < T:
    while t < T_E:
        # 2.3更新劣化水平、故障率和库存等
        ε = np.random.uniform(0, 0.5, 1)  # 采样时间
        X_increment = np.random.gamma(α_X * ε, β_X)  # 随机劣化增量
        X += X_increment
        t += ε
        # 更新故障率和采样分数等
        base0 = λ * X
        if base0 == 0:
            # 处理底部为零的情况，给 base 加上一个小值来避免除以零
            base0 += epsilon0
        base1 = α * θ * (X ** (θ - 1)) + β * (γ + λ * X) * (X ** (γ - 1)) * np.expm1(base0)
        if base1 == 0:
            # 处理底部为零的情况，给 base 加上一个小值来避免除以零
            base1 += epsilon1
        p = p_0 + p_1 * (1 - np.expm1(-base1))
        f = 1 - np.exp(- Is * base1)
        m = m_0 + m_1 * (1 - np.expm1(-base1))
        # 更新产品数量和库存
        if n_oldsystem == 0:
            if I >= Z:
                N_def += (D * f * p * m * ε) / (1 - f * p * m)
                N_rep += (D * f * p * ε) / (1 - f * p * m)
                N_ins += (D * f * (1 + p) * ε) / (1 - f * p * m)
                C_h += c_h * Z * ε
                N_tcb += (D * p * (1 - f) * ε) / (1 - f * p * m)
                N_tc += D * ε
                I = Z
            else:
                N_def += R_max * f * p * m * ε
                N_rep += R_max * f * p * ε
                N_ins += R_max * f * (1 + p) * ε
                C_h += c_h * (R_max * (1 - f * p * m) - D) * ε
                N_tcb += R_max * p * (1 - f) * ε
                N_tc += D * ε
                I += (R_max * (1 - f * p * m) - D) * ε
        else:
            if I >= Z:
                N_def += (D * f * p * m * ε) / (1 - f * p * m)
                N_rep += (D * f * p * ε) / (1 - f * p * m)
                N_ins += (D * f * (1 + p) * ε) / (1 - f * p * m)
                C_h += c_h * Z * ε
                N_tcb += (D * p * (1 - f) * ε) / (1 - f * p * m)
                N_tc += D * ε
                N_tcb_new += (D * p * (1 - f) * ε) / (1 - f * p * m)
                I = Z
            else:
                N_def += R_max * f * p * m * ε
                N_rep += R_max * f * p * ε
                N_ins += R_max * f * (1 + p) * ε
                C_h += c_h * (R_max * (1 - f * p * m) - D) * ε
                N_tcb += R_max * p * (1 - f) * ε
                N_tc += D * ε
                N_tcb_new += R_max * p * (1 - f) * ε
                I += (R_max * (1 - f * p * m) - D) * ε
        # 2.4将劣化水平X与故障阈值L进行比较
        if X >= L:
            if I >= Z:
                # 3.1
                # 实施CM并计算第n个批次产生的总成本
                # Scenario3
                n_3 += 1
                I = Z
                T_CM = np.random.exponential(scale=1.0 / 1.21, size=1)  # CM持续时间
                C_set = c_set  # 生产运行的设置成本
                C_h += c_h * (Z ** 2) / (2 * D)
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                if T_CM >= (Z / D):
                    C_s = c_s * 0.5 * D * ((T_CM - Z / D) ** 2)
                    k += 1
                else:
                    C_s = 0
                C_Total_3 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_CM + C_s
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                CM += c_CM
                s += C_s
                C_Total += C_Total_3
            else:
                # 3.2
                # 实施CM并计算第n个批次产生的总成本
                # Scenario1
                n_1 += 1
                T_CM = np.random.exponential(scale=1.0 / 1.21, size=1)  # CM持续时间
                C_set = c_set  # 生产运行的设置成本
                C_h += c_h * (I ** 2) / (2 * D)
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                if T_CM >= (I / D):
                    C_s = c_s * 0.5 * D * ((T_CM - I / D) ** 2)
                    k += 1
                else:
                    C_s = 0
                C_Total_1 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_CM + C_s
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                CM += c_CM
                s += C_s
                C_Total += C_Total_1
            # 3.3
            # 更新劣化水平、Gamma分布参数、库存水平和总时间
            if T_CM >= (I / D):
                t += T_CM
            else:
                t += I / D
            I = 0
            T_Total += t
            t = 0
            X = 0
            N_def = 0  # 当前生产的缺陷产品数
            N_rep = 0  # 当前生产的需要返工的产品数
            N_ins = 0  # 到目前为止的检查次数
            N_tcb = 0  # 到目前为止发送给顾客的坏品数
            N_tc = 0  # 到目前为止发送给顾客的产品数
            n_oldsystem = 0
            AOQ = 0
            C_h = 0
            if t == 0:
                break
        else:
            # 2.5
            # 计算AOQ
            # AOQ与阈值AOQL进行比较
            AOQ = (N_tcb / N_tc)
            if AOQ >= AOQL:
                if I >= Z:
                    # 4.1
                    # 实施大修并计算第n个批次产生的总成本
                    # Scenario4
                    n_4 += 1
                    I = Z
                    T_OH = np.random.exponential(scale=1.0 / 0.95, size=1)  # OH持续时间
                    C_set = c_set  # 生产运行的设置成本
                    C_h += c_h * (Z ** 2) / (2 * D)
                    C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                    C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                    C_ins = c_ins * N_ins  # 单位检查成本
                    if n_oldsystem == 0:
                        C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                    else:
                        C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                    if T_OH >= (Z / D):
                        C_s = c_s * 0.5 * D * ((T_OH - Z / D) ** 2)
                        k += 1
                    else:
                        C_s = 0
                    C_Total_4 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_OH + C_s
                    set += C_set
                    h += C_h
                    rej += C_rej
                    rep += C_rep
                    ins += C_ins
                    q += C_q
                    OH += c_OH
                    s += C_s
                    C_Total += C_Total_4
                else:
                    # 4.2
                    # 实施大修并计算第n个批次产生的总成本
                    # Scenario2
                    n_2 += 1
                    T_OH = np.random.exponential(scale=1.0 / 0.95, size=1)  # OH持续时间
                    C_set = c_set  # 生产运行的设置成本
                    C_h += c_h * (I ** 2) / (2 * D)
                    C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                    C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                    C_ins = c_ins * N_ins  # 单位检查成本
                    if n_oldsystem == 0:
                        C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                    else:
                        C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                    if T_OH >= (I / D):
                        C_s = c_s * 0.5 * D * ((T_OH - I / D) ** 2)
                        k += 1
                    else:
                        C_s = 0
                    C_Total_2 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_OH + C_s
                    set += C_set
                    h += C_h
                    rej += C_rej
                    rep += C_rep
                    ins += C_ins
                    q += C_q
                    OH += c_OH
                    s += C_s
                    C_Total += C_Total_2
                # 4.3
                # 更新劣化水平、Gamma分布参数、库存水平和总时间
                if T_OH >= (Z / D):
                    t += T_OH
                else:
                    t += I / D
                I = 0
                T_Total += t
                t = 0
                X = 0
                N_def = 0  # 当前生产的缺陷产品数
                N_rep = 0  # 当前生产的需要返工的产品数
                N_ins = 0  # 到目前为止的检查次数
                N_tcb = 0  # 到目前为止发送给顾客的坏品数
                N_tc = 0  # 到目前为止发送给顾客的产品数
                n_oldsystem = 0
                AOQ = 0
                C_h = 0
                if t == 0:
                    break
    while t >= T_E:
        if X >= L:
            if I >= Z:
                # 3.1
                # 实施CM并计算第n个批次产生的总成本
                # Scenario3
                print('!!!!!!')
                n_3 += 1
                I = Z
                T_CM = np.random.exponential(scale=1.0 / 1.21, size=1)  # CM持续时间
                C_set = c_set  # 生产运行的设置成本
                C_h += c_h * (Z ** 2) / (2 * D)
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                if T_CM >= (Z / D):
                    C_s = c_s * 0.5 * D * ((T_CM - Z / D) ** 2)
                    k += 1
                else:
                    C_s = 0
                C_Total_3 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_CM + C_s
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                CM += c_CM
                s += C_s
                C_Total += C_Total_3
            else:
                # 3.2
                # 实施CM并计算第n个批次产生的总成本
                # Scenario1
                n_1 += 1
                T_CM = np.random.exponential(scale=1.0 / 1.21, size=1)  # CM持续时间
                C_set = c_set  # 生产运行的设置成本
                C_h += c_h * (I ** 2) / (2 * D)
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                if T_CM >= (I / D):
                    C_s = c_s * 0.5 * D * ((T_CM - I / D) ** 2)
                    k += 1
                else:
                    C_s = 0
                C_Total_1 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_CM + C_s
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                CM += c_CM
                s += C_s
                C_Total += C_Total_1
            # 3.3
            # 更新劣化水平、Gamma分布参数、库存水平和总时间
            if T_CM >= (I / D):
                t += T_CM
            else:
                t += I / D
            I = 0
            T_Total += t
            t = 0
            X = 0
            N_def = 0  # 当前生产的缺陷产品数
            N_rep = 0  # 当前生产的需要返工的产品数
            N_ins = 0  # 到目前为止的检查次数
            N_tcb = 0  # 到目前为止发送给顾客的坏品数
            N_tc = 0  # 到目前为止发送给顾客的产品数
            n_oldsystem = 0
            AOQ = 0
            C_h = 0
            if t == 0:
                break
        # 2.5
        # 计算AOQ
        # AOQ与阈值AOQL进行比较
        elif (N_tcb / N_tc) >= AOQL:
            if I >= Z:
                # 4.1
                # 实施大修并计算第n个批次产生的总成本
                # Scenario4
                n_4 += 1
                I = Z
                T_OH = np.random.exponential(scale=1.0 / 0.95, size=1)  # OH持续时间
                C_set = c_set  # 生产运行的设置成本
                C_h += c_h * (Z ** 2) / (2 * D)
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                if T_OH >= (Z / D):
                    C_s = c_s * 0.5 * D * ((T_OH - Z / D) ** 2)
                    k += 1
                else:
                    C_s = 0
                C_Total_4 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_OH + C_s
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                OH += c_OH
                s += C_s
                C_Total += C_Total_4
            else:
                # 4.2
                # 实施大修并计算第n个批次产生的总成本
                # Scenario2
                n_2 += 1
                T_OH = np.random.exponential(scale=1.0 / 0.95, size=1)  # OH持续时间
                C_set = c_set  # 生产运行的设置成本
                C_h += c_h * (I ** 2) / (2 * D)
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                if T_OH >= (I / D):
                    C_s = c_s * 0.5 * D * ((T_OH - I / D) ** 2)
                    k += 1
                else:
                    C_s = 0
                C_Total_2 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_OH + C_s
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                OH += c_OH
                s += C_s
                C_Total += C_Total_2
            # 4.3
            # 更新劣化水平、Gamma分布参数、库存水平和总时间
            if T_OH >= (Z / D):
                t += T_OH
            else:
                t += I / D
            I = 0
            T_Total += t
            t = 0
            X = 0
            N_def = 0  # 当前生产的缺陷产品数
            N_rep = 0  # 当前生产的需要返工的产品数
            N_ins = 0  # 到目前为止的检查次数
            N_tcb = 0  # 到目前为止发送给顾客的坏品数
            N_tc = 0  # 到目前为止发送给顾客的产品数
            n_oldsystem = 0
            AOQ = 0
            C_h = 0
            if t == 0:
                break
        else:
            # 步骤5
            if X >= ξ:
                # 实施PM并计算第n个批次产生的总成本
                # Scenario5
                n_5 += 1
                # I = Z
                T_PM = np.random.exponential(scale=1.0 / 0.79, size=1)  # PM持续时间
                C_set = c_set  # 生产运行的设置成本
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                # 计算总成本
                # 更新劣化水平、Gamma分布参数、库存水平和总时间
                X_new = beta.rvs(u, v) * ξ
                if T_PM <= I / D:
                    C_h += c_h * (I ** 2) / (2 * D)
                    C_s = 0
                    t += I / D
                    I = 0
                else:
                    C_h += c_h * (I ** 2) / (2 * D)
                    C_s = 0.5 * D * ((T_PM - I / D) ** 2)
                    k += 1
                    t += T_PM
                    I = 0
                C_Total_5 = C_set + C_h + C_rej + C_rep + C_ins + C_q + c_PM + C_s
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                PM += c_PM
                s += C_s
                C_Total += C_Total_5
                T_Total += t
                t = 0
                X = X_new
                α_X *= g
                N_def = 0  # 当前生产的缺陷产品数
                N_rep = 0  # 当前生产的需要返工的产品数
                N_ins = 0  # 到目前为止的检查次数
                N_tcb = 0  # 到目前为止发送给顾客的坏品数
                N_tc = 0  # 到目前为止发送给顾客的产品数
                n_oldsystem = 0
                AOQ = 0
                C_h = 0
                if t == 0:
                    break
            else:
                # 不进行任何操作，结束循环并计算第n个批次产生的总成本
                # Scenario6
                n_6 += 1
                # I = Z
                C_h += c_h * (I ** 2) / (2 * D)
                C_set = c_set  # 生产运行的设置成本
                C_rej = c_rej * N_def  # 返工缺陷产品的拒绝成本
                C_rep = c_rep * N_rep  # 缺陷产品的维修成本
                C_ins = c_ins * N_ins  # 单位检查成本
                if n_oldsystem == 0:
                    C_q = c_q * N_tcb  # 出售给客户不良商品的单位成本
                else:
                    C_q = c_q * N_tcb_new  # 出售给客户不良商品的单位成本
                C_Total_6 = C_set + C_h + C_rej + C_rep + C_ins + C_q
                set += C_set
                h += C_h
                rej += C_rej
                rep += C_rep
                ins += C_ins
                q += C_q
                C_Total += C_Total_6
                t += I / D
                # 更新库存水平和总时间
                I = 0
                T_Total += t
                t = 0
                N_def = 0  # 当前生产的缺陷产品数
                N_rep = 0  # 当前生产的需要返工的产品数
                N_ins = 0  # 到目前为止的检查次数
                N_tcb_new = 0  # 到目前为止发送给顾客的坏品数
                C_h = 0
                n_oldsystem += 1
                if t == 0:
                    break
# 步骤6
# 计算并输出ETC
ETC = C_Total / T_Total
print(Z, ' ', T_E, ' ', ξ, ' ', Is, ' ', ETC)

'''print('安全库存时刻在循环结束之后的次数为：', num)
print('ETC=', ETC)
print('C_Total=', C_Total)
print('T_Total=', T_Total)
print('n_1=', n_1)
print('n_2=', n_2)
print('n_3=', n_3)
print('n_4=', n_4)
print('n_5=', n_5)
print('n_6=', n_6)
print(ETC, C_Total, T_Total, n_1, n_2, n_3, n_4, n_5, n_6)

print(a, b, c, d, e, f)
print(values_T_E, values_Z, values_S, values_AOQL, values_ξ, values_Is)'''

print('此时ETC为：', ETC, '场景1-6出现的次数分别为：', n_1, ',', n_2, ',', n_3, ',', n_4, ',', n_5, ',', n_6)

'''end_time = time.time()  # 获取当前时间

print(f"The code ran in {end_time - start_time} seconds")
print()'''