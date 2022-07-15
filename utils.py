import os
import re
import math
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def distance(v1, v2):
    '''
    :param v1, v2:要求都是向量(numpy)
    :return: 两点的距离
    '''
    return np.linalg.norm(v1 - v2)

def radialDistributionFunction(pos, step = 0.1):
    '''
    :param pos:三维坐标(x,3)(numpy)
    :return:RDF分布
    '''
    N = pos.shape[0]
    center = np.average(pos, axis = 0)
    print(f'中心点为{center}')
    dis = np.array([])
    for i in range(N):
        dis = np.append(dis, distance(pos[i], center))
    maxDistance = np.max(dis)
    density = np.divide(N, 4 / 3 * np.pi * np.float_power(maxDistance, 3)) # 采用球形体积计算密度
    x = np.arange(start = 0, stop = maxDistance, step = step)
    y = np.zeros_like(x)
    for i in range(N):
        y[int(dis[i] / step)] += 1
    # print(y)
    for i in range(N):
        s_i = 2 * np.pi * x[i] * step
        if s_i > 0:
            y[i] = np.divide(y[i], s_i * N * density)
        else:
            y[i] = 0
    # print(y)
    return x, y