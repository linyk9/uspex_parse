import os
import re
import math
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

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
    return x, y

def POSCARSplit(filename, output = r'.\POSCAR'):
    '''
    分割uspex运行出的总POSCAR文件
    :param filename: 一般是***POSCAR文件
    :param output: 输出目录
    :return:
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[:2] == 'EA':
                try:
                    fTmp.close()
                except:
                    pass
                name = line.replace('\n', '')
                tmpPath = os.path.join(output, rf'{name}')
                mkdir(tmpPath)
                if not os.path.isfile(rf'{tmpPath}\POSCAR'):
                    fTmp = open(rf'{tmpPath}\POSCAR', 'w')
            try:
                fTmp.write(line)
            except:
                pass
        try:
            fTmp.close()
        except:
            pass
    print('Finished!')


def readPOSCAR(filename):
    '''
    读入POSCAR文件
    :param filename: POSCAR文件名
    :return: Angstrom坐标数组(numpy)
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
        latticeVectors = []
        atoms = []
        for index, line in enumerate(lines):
            if index in [0, 1, 5, 6, 7]:
                # 多余行处理
                continue
            elif index in [2, 3, 4]:
                # 晶格参数
                line = line.replace('\n', '').split(' ')
                line = [l for l in line if l]
                latticeVectors.append(line)
            else:
                # 晶格坐标
                line = line.replace('\n', '').split(' ')
                line = [l for l in line if l]
                atoms.append(line)
        latticeVectors = np.array(latticeVectors, dtype = np.float)
        atoms = np.array(atoms, dtype = np.float)
    # print(atoms)
    # print(latticeVectors)
    return np.matmul(atoms, latticeVectors)

def drawRDF(filename):
    pos = readPOSCAR(r'.\POSCAR\EA114\POSCAR')
    print(pos)
    r, gr = radialDistributionFunction(pos)
    plt.plot(r, gr, lw = '1', c = 'r')
    plt.title('radial distribution function')
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.show()