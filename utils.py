import os
import re
import math
from decimal import Decimal
from functools import singledispatch
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

def num_density(pos, maxDistance):
    '''
    采用球形体积计算密度
    '''
    N = pos.shape[0]
    density = np.divide(N, 4 / 3 * np.pi * np.float_power(maxDistance, 3)) # 采用球形体积计算密度
    return density

def radialDistributionFunction(pos, step = 0.01, maxDistance = 10):
    '''
    :param pos:三维坐标(x,3)(numpy)
    :return:RDF分布
    '''
    N = pos.shape[0]
    #center = np.average(pos, axis = 0)
    center = pos[0]
    print(f'中心点为{center}')
    dis = np.array([])
    for i in range(N):
        dis = np.append(dis, distance(pos[i], center))
    print(f'每个原子的距离{dis}')
    density = num_density(pos,maxDistance)
    x = np.arange(start = 0, stop = maxDistance, step = step)
    gr = np.zeros_like(x)

    # 看看每个原子落在哪一个小区间上
    for i in range(N):
        gr[int(dis[i] / step)] += 1

    y = []
    for i in range(len(x)):
        s_i = 4 * np.pi * x[i] * x[i]
        # rho : 局部密度
        if s_i > 0:
            rho = np.divide(gr[i], s_i)
        else:
            rho = 0
        y.append(rho / density)

    # integral (dens * g(r) * 4*pi*r^2) ==> N - 1
    sum_y = 0
    for i in range(len(y)):
        sum_y += y[i] * 4 * np.pi * x[i] * x[i] * density
    print(f'integral={sum_y}(N - 1)')
    return x, y 


def POSCARSplit(filename, output = r'./POSCAR'):
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
    return np.matmul(atoms, latticeVectors)

def drawRDF(type = 0, arg = [], isPrint = False):
    '''
    :param type: 0代表传入POSCAR文件计算；1代表直接给r和gr计算
    :param arg: 列表变量；0表示传入文件名；1表示传入r，gr，保存的文件名
    :param isPrint: 是否展示可视化图像
    '''
    if type == 0:
        filename = arg[0]
        pos = readPOSCAR(filename)
        maxdis = 10
        step = 0.01
        r, gr = radialDistributionFunction(pos, step, maxdis)
    else:
        r, gr, filename = arg[0], arg[1], arg[2]
    plt.plot(r, gr, lw = '1', c = 'r')
    plt.title('radial distribution function')
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.savefig(f'{filename}.png')
    if isPrint:
        plt.show()

def drawRDFs(dir, isPrint = False):
    for d in os.listdir(dir):
        try:
            print(f'当前处理文件：{d}')
            drawRDF(0, [f'{dir}\\{d}'], isPrint)
            print()
        except:
            print(f'{d}文件处理失败\n')