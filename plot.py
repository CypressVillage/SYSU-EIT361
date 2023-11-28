import os
import matplotlib.pyplot as plt

actions = (
    # 'calculate',
    'plot',
)

decode_methods = (
    'viterbi_hard',
    'viterbi_soft',
    'bcjr',
    'turbo',
)

# 判断assets文件夹是否存在，不存在则创建
if not os.path.exists('assets'):
    os.mkdir('assets')

if 'calculate' in actions:
    for method in decode_methods:
        # 如果存在data_method.txt文件，则删除
        if os.path.exists('assets/data_' + method + '.txt'):
            os.remove('assets/data_' + method + '.txt')

        # 在assets文件夹中运行main.exe
        if os.name == 'posix':
            os.system('cd assets && ../build/main ' + method)
        else:
            os.system('cd assets && ..\\build\main.exe ' + method)

        os.rename('assets/data.txt', 'assets/data_' + method + '.txt')

if 'plot' in actions:
    plt.figure()
    for method in decode_methods:
        snr = []
        ber = []
        with open('assets/data_'+method+'.txt', 'r') as f:
            for line in f.readlines():
                snr.append(float(line.split()[0]))
                ber.append(float(line.split()[1]))
        plt.semilogy(snr, ber, '-o', label=method)
    plt.xlabel('SNR(dB)')
    plt.ylabel('BER')
    plt.grid()
    plt.legend()
    plt.savefig('assets/figure.png')
