import os
import matplotlib.pyplot as plt

# 以不同的参数运行build/main.exe
decode_methods = [
    'viterbi_hard',
    'viterbi_soft',
    'bcjr',
    'turbo',
]

# 判断assets文件夹是否存在，不存在则创建
if not os.path.exists('assets'):
    os.mkdir('assets')

for method in decode_methods:
    os.system('build\main.exe ' + method)

    # 读取data.txt文件中的SNR和BER数据
    snr = []
    ber = []
    with open('assets/data.txt', 'r') as f:
        for line in f.readlines():
            snr.append(float(line.split()[0]))
            ber.append(float(line.split()[1]))
    os.rename('assets/data.txt', 'assets/data_' + method + '.txt')

    # 绘制SNR-BER曲线
    plt.figure(figsize=(8, 6))
    plt.semilogy(snr, ber, 'b-o')
    plt.xlabel('SNR(dB)')
    plt.ylabel('BER')
    plt.grid(True)
    plt.title(method)
    plt.savefig('assets/figure_' + method + '.png')