# 5个参数：α, β, θ, γ, λ
import numpy as np
import matplotlib.pyplot as plt
# 坐标轴负号
import matplotlib
matplotlib.rcParams['axes.unicode_minus'] = False

α = 1  # 浴盆曲线的尺度参数
β = 0.00025  # 浴盆曲线的尺度参数
θ = 0.05  # 浴盆曲线的形状参数
γ = 0.01  # 浴盆曲线的形状参数
λ = 1.25  # 浴盆曲线的加速度参数
p_0 = 0.02  # 缺陷率的初值
p_1 = 0.28  # 缺陷率的峰值

x = np.linspace(0, 10, 1000)
p = p_0 + p_1 * (1 - np.exp(- (α * θ * (x ** (θ - 1)) + β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x))))
plt.figure(figsize=(10, 4))
# 中文标签
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.plot(x, p, color='blue', linewidth=1, linestyle='-')
plt.xlabel('X(t)')
plt.ylabel('Quality level p(X)')
plt.title('Quality degradation with changes in degradation level X(t)')
plt.grid(False)
plt.show()
