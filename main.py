import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# print(__version__)
def main():
    # 1.获取好的结构,并写入csv
    filename = './uspex/results1/goodStructures'
    with open(filename, 'r') as f:
        text = f.read()
    text = text.replace(' ', '')
    # 正则表达式匹配出每一行数据
    pattern = re.compile(r"\|(\d+)\|\d+\|(\w+)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|\n")
    match = pattern.findall(string = text)
    print(f'共有 {len(match)} 个较好结构')
    columns = ['ID', 'Origin', 'Composition', 'Enthalpy (eV)', 'Volume (A^3)', 'SYMMETRY (N)', 'Structure order', 'Average order', 'Quasientropy']
    db = pd.DataFrame(data = match, columns = columns)
    print(db.head())
    db.to_csv('.\goodStructures.csv')

    # 2.获取较稳定和亚稳定结构
    stable = 0.025852 # 1 / 2 k_b * t 在300k下为0.025852 eV
    metastable = 0.05
    db2 = pd.DataFrame()
    db3 = pd.DataFrame()
    bestEnthalpy = float(db['Enthalpy (eV)'][0])
    for index, value in enumerate(db['Enthalpy (eV)']):
        value = float(value)
        if value <= bestEnthalpy + stable:
            db2 = db2.append(db[db.index == index])

        if value <= bestEnthalpy + metastable:
            db3 = db3.append(db[db.index == index])

    db2.to_csv('.\stableStructures.csv')
    db3.to_csv('.\metastableStructures.csv')



if __name__ == '__main__':
    main()