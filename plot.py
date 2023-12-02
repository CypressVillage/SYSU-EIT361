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
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    plt.figure()
    msg_len, repeat_times, iter_times = 0, 0, 0

    for method in decode_methods:
        snr = []
        ber = []
        with open('assets/data_'+method+'.txt', 'r') as f:
            # 第一行是message_length，重复次数，turbo译码迭代次数
            msg_len, repeat_times, iter_times = f.readline().split()
            msg_len, repeat_times, iter_times = int(msg_len), int(repeat_times), int(iter_times)

            for line in f.readlines():
                snr.append(float(line.split()[0]))
                ber.append(float(line.split()[1]))
        plt.semilogy(snr, ber, '-o', label=method)

    infostr = f'''
    信息序列长度: {msg_len}
    重复次数: {repeat_times}
    turbo译码迭代次数: {iter_times}
    '''
    plt.text(0.6, 0.5, infostr, fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel('SNR(dB)')
    plt.ylabel('BER')
    plt.grid()
    plt.legend()
    plt.savefig('assets/figure.png')
