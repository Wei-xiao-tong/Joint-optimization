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
Is1 = 1000
Is2 = 100
Is3 = 50
Is4 = 10
Is5 = 5
Is6 = 2
Is7 = 0.5
Is8 = 0

x = np.linspace(0, 10, 1000)
p = p_0 + p_1 * (1 - np.exp(- (α * θ * (x ** (θ - 1)) + β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x))))
f1 = 1 - np.exp(- (α * θ * Is1 * (x ** (θ - 1)) + Is1 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
f2 = 1 - np.exp(- (α * θ * Is2 * (x ** (θ - 1)) + Is2 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
f3 = 1 - np.exp(- (α * θ * Is3 * (x ** (θ - 1)) + Is3 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
f4 = 1 - np.exp(- (α * θ * Is4 * (x ** (θ - 1)) + Is4 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
f5 = 1 - np.exp(- (α * θ * Is5 * (x ** (θ - 1)) + Is5 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
f6 = 1 - np.exp(- (α * θ * Is6 * (x ** (θ - 1)) + Is6 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
f7 = 1 - np.exp(- (α * θ * Is7 * (x ** (θ - 1)) + Is7 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
f8 = 1 - np.exp(- (α * θ * Is8 * (x ** (θ - 1)) + Is8 * β * (γ + λ * x) * (x ** (γ - 1)) * np.exp(λ * x)))
plt.figure(figsize=(10, 4))
# 中文标签
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.plot(x, p, color='red', linewidth=1, linestyle='-', label='p(X)')
plt.plot(x, f1, color='orange', linewidth=1, linestyle='-', label='Is=1000')
plt.plot(x, f2, color='orange', linewidth=1, linestyle='-.', label='Is=100')
plt.plot(x, f3, color='green', linewidth=1, linestyle='-', label='Is=50')
plt.plot(x, f4, color='green', linewidth=1, linestyle='-.', label='Is=10')
plt.plot(x, f5, color='blue', linewidth=1, linestyle='-', label='Is=5')
plt.plot(x, f6, color='blue', linewidth=1, linestyle='-.', label='Is=2')
plt.plot(x, f7, color='black', linewidth=1, linestyle='-', label='Is=0.5')
plt.plot(x, f8, color='black', linewidth=1, linestyle='-.', label='Is=0')
plt.xlabel('X(t)')
plt.ylabel('inspection fraction')
plt.title('The variation of the inspected fraction with different values of Is')
plt.grid(False)
plt.legend(loc="upper right")
plt.show()
