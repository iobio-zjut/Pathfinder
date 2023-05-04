"""
!/usr/bin/env python
 -*- coding: utf-8 -*-
 @Time    : 2022/3/30 14:03
 @Author  : hzh
 @File    : ContactOreder.py
"""
import argparse
import csv
import math
import os

# from utils import Contact_Order, read_pdb
import sys


def read_pdb(f_path):
    """
        read .pdb file get CA,N,C to rebuild CB coordination
        CB is used to calculate contact

    :param f_path: path of pdb
    :return:
    """
    pdb_N = []
    pdb_CA = []
    pdb_C = []
    with open(f_path, 'r') as f1:
        reader = csv.reader(f1)
        for row in reader:
            row = row[0]
            if row.split()[0] == 'TER':
                break

            if row.split()[2] == 'N':
                temp = [float(row[31:38]),
                        float(row[39:46]),
                        float(row[47:54])]
                pdb_N.append(temp)

            if row.split()[2] == 'CA':
                temp = [float(row[31:38]),
                        float(row[39:46]),
                        float(row[47:54])]
                pdb_CA.append(temp)

            if row.split()[2] == 'C':
                temp = [float(row[31:38]),
                        float(row[39:46]),
                        float(row[47:54])]
                pdb_C.append(temp)

    return pdb_N, pdb_CA, pdb_C


def Contact_Order(pdb_path):
    """
        Contact order equation : (1 / N*L) * sum(sequence separation)
        sequence separation is the index of i and j distance in sequence

        In original paper contact use CB,and use CA insteal of CB,when
        there is no CB. So there is a little different from the data of
        the paper.

        'Plaxco K W, Simons K T, Baker D. Contact order,
        transition state placement and the refolding rates of single domain proteins[J].
        Journal of molecular biology, 1998, 277(4): 985-994.'

    :param pdb_N:
    :param pdb_CA:
    :param pdb_C:
    :return: CO: contact order
    """

    pdb_N, pdb_CA, pdb_C = read_pdb(pdb_path)

    ca = -0.58273431
    cb = 0.56802827
    cc = -0.54067466

    L = len(pdb_CA)
    deltaS = 0
    N = 0
    pdb_CB = []
    res_CO = []

    # rebuild CB
    for i in range(L):
        b = [pdb_CA[i][0] - pdb_N[i][0],
             pdb_CA[i][1] - pdb_N[i][1],
             pdb_CA[i][2] - pdb_N[i][2]]
        c = [pdb_C[i][0] - pdb_CA[i][0],
             pdb_C[i][1] - pdb_CA[i][1],
             pdb_C[i][2] - pdb_CA[i][2]]
        a = [b[1] * c[2] - b[2] * c[1],
             b[2] * c[0] - b[0] * c[2],
             b[0] * c[1] - b[1] * c[0]]
        temp = [ca * a[0] + cb * b[0] + cc * c[0] + pdb_CA[i][0],
                ca * a[1] + cb * b[1] + cc * c[1] + pdb_CA[i][1],
                ca * a[2] + cb * b[2] + cc * c[2] + pdb_CA[i][2]]
        pdb_CB.append(temp)

    for i in range(L):
        temp = 0
        for j in range(L):
            if i == j:
                continue
            x = pdb_CB[i][0] - pdb_CB[j][0]
            y = pdb_CB[i][1] - pdb_CB[j][1]
            z = pdb_CB[i][2] - pdb_CB[j][2]

            dist = math.sqrt(x * x + y * y + z * z)

            # if < 8A , contact
            if dist < 8:
                N += 1
                deltaS += abs(i - j)

            if dist < 20:
                temp += abs(i - j) / dist

        res_CO.append(temp)

    return deltaS / (L * N * 1.0), res_CO


if __name__ == '__main__':
    # add argument
    parser = argparse.ArgumentParser(description="calculate contact order.")
    parser.add_argument('-p', type=str, default='./', help='path of working directory.')
    args = parser.parse_args()

    path = args.p
    file = os.listdir(path)

    file = [i for i in file if i[-4:] == '.pdb']

    if 'native.pdb' not in file:
        sys.exit("lack native")
    file.remove('native.pdb')
    NCO, Nresi = Contact_Order(path + '/native.pdb')
    CO = {'native  ': round(NCO * 100, 1)}

    for file_index in file:

        t_CO, resi = Contact_Order(path + "/" + file_index)
        CO[file_index[:-4]] = (round(t_CO * 100, 1))
        result = []
        CA = 0
        with open(path + "/" + file_index, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                result.append(row[0])
                if row[0].split()[0] == 'TER':
                    break
        for index in range(len(result)):

            if result[index].split()[0] == 'TER':
                break
            replace_resi = round((1 - min(resi[CA] / Nresi[CA], 1)) * 100, 1)
            if 0 < replace_resi < 10:
                replace_resi = "  " + str(replace_resi)
            elif replace_resi == 0:
                replace_resi = '  0.00'
            else:
                replace_resi = " " + str(replace_resi)

            result[index] = result[index][:61] + replace_resi + result[index][66:]
            if result[index].split()[2] == 'CA':
                if CA < len(resi) - 1:
                    CA += 1

        with open(path + "/" + file_index, 'w') as w:
            for i in result:
                w.write(i + '\n')

    with open(path + '/ContactOrder.txt', 'w') as w:
        rCO = [i for i in CO.items()]
        rCO = sorted(rCO, key=lambda x: x[1], reverse=True)

        print(rCO)

        for index in rCO:
            w.write(index[0] + '\t\t' + str(index[1]) + '\n')
