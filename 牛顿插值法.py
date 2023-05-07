# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 14:42:25 2023

@author: 86319
"""

# 定义温度和表面张力系数列表
T = [0, 5, 10, 15, 20, 25, 30]
mu = [75.64, 74.92, 74.22, 73.49, 72.75, 71.97, 71.18]

# 定义需要插值的温度值
T_interp = 13.2

# 三点拉格朗日插值法
def lagrange_interp(x, x_list, y_list):
    n = len(x_list)
    assert n == len(y_list)
    assert n >= 3
    if x <= x_list[1]:
        i = 0
    elif x >= x_list[n-2]:
        i = n-3
    else:
        for j in range(1, n-2):
            if x >= x_list[j] and x <= x_list[j+1]:
                i = j
                break
    x0, x1, x2 = x_list[i:i+3]
    y0, y1, y2 = y_list[i:i+3]
    L0 = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))
    L1 = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
    L2 = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1))
    y_interp = y0 * L0 + y1 * L1 + y2 * L2
    return y_interp

# 计算13.2℃处的表面张力系数的近似值
mu_interp_lagrange = lagrange_interp(T_interp, T, mu)

# 输出结果
print("三点拉格朗日插值法：")
print("13.2℃时的表面张力系数的近似值为：{:.8f}".format(mu_interp_lagrange))

# 牛顿插值法
def newton_interpolation(x, T=T, mu=mu):
    n = len(T) - 1
    diff_quot = [mu]
    for i in range(n):
        diff_quot_i = []
        for j in range(n - i):
            diff_quot_i.append((diff_quot[i][j + 1] - diff_quot[i][j]) / (T[j + i + 1] - T[j]))
        diff_quot.append(diff_quot_i)
    p = diff_quot[0][0]
    for i in range(1, n + 1):
        term = diff_quot[i][0]
        for j in range(i):
            term *= (x - T[j])
        p += term
    return p

# 计算13.2℃时的表面张力系数
x = 13.2
p = newton_interpolation(x)

# 输出结果
print("牛顿插值法：")
print("13.2℃时的表面张力系数的近似值为：{:.8f}".format(p))
