import os
import matplotlib.pyplot as plt

actions = (
    # 'calculate',      # 计算并保存数据
    # 'plot',           # 画图
    'calculate_turbo',  # 计算并保存不同迭代次数turbo译码数据
    'plot_turbo'        # 画不同迭代次数turbo译码图
)

decode_methods = ( # 仅calculate有效
    'viterbi_hard',
    'viterbi_soft',
    'bcjr',
    'turbo',
)

turbo_iter_times = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] # 仅calculate_turbo有效

is_linux = os.name == 'posix'

# 判断assets文件夹是否存在，不存在则创建
if not os.path.exists('assets'):
    os.mkdir('assets')

if 'calculate' in actions:
    for method in decode_methods:
        # 如果存在data_method.txt文件，则删除
        if os.path.exists('assets/data_' + method + '.txt'):
            os.remove('assets/data_' + method + '.txt')

        # 在assets文件夹中运行main.exe
        if is_linux:
            os.system('cd assets && ../build/main ' + method)
        else:
            os.system('cd assets && ..\\build\main.exe ' + method)

        os.rename('assets/data.txt', 'assets/data_' + method + '.txt')

if 'plot' in actions:
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    plt.figure()
    msg_len, repeat_times, iter_times, conv_1, conv_2 = 0, 0, 0, 0, 0

    for method in decode_methods:
        snr = []
        ber = []
        with open('assets/data_'+method+'.txt', 'r') as f:
            # 第一行是message_length，重复次数，turbo译码迭代次数
            msg_len, repeat_times, iter_times, conv_1, conv_2 = f.readline().split()
            msg_len, repeat_times, iter_times, conv_1, conv_2 = int(msg_len), int(repeat_times), int(iter_times), int(conv_1), int(conv_2)
            # 从第二行开始是snr和ber
            for line in f.readlines():
                snr.append(float(line.split()[0]))
                ber.append(float(line.split()[1]))
        plt.semilogy(snr, ber, '-o', label=method)

    infostr = f'''
    信息序列长度: {msg_len}
    重复次数: {repeat_times}
    turbo译码迭代次数: {iter_times}
    卷积码类型: ({conv_1}, {conv_2})
    '''
    plt.text(0.6, 0.5, infostr, fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel('SNR(dB)')
    plt.ylabel('BER')
    plt.grid()
    plt.legend()
    plt.savefig('assets/figure.png')

if 'calculate_turbo' in actions:
    for iter_time in turbo_iter_times:
        # 如果存在data_turbo_iter_time.txt文件，则删除
        if os.path.exists('assets/data_turbo_iter_'+str(iter_time)+'.txt'):
            os.remove('assets/data_turbo_iter_'+str(iter_time)+'.txt')

        # 在assets文件夹中运行main.exe
        if is_linux:
            os.system('cd assets && ../build/main turbo ' + str(iter_time))
        else:
            os.system('cd assets && ..\\build\main.exe turbo ' + str(iter_time))

        os.rename('assets/data.txt', 'assets/data_turbo_iter_'+str(iter_time)+'.txt')

if 'plot_turbo' in actions:
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    plt.figure()
    msg_len, repeat_times, iter_times, conv_1, conv_2 = 0, 0, 0, 0, 0

    for iter_time in turbo_iter_times:
        snr = []
        ber = []
        with open('assets/data_turbo_iter_'+str(iter_time)+'.txt', 'r') as f:
            # 第一行是message_length，重复次数，turbo译码迭代次数
            msg_len, repeat_times, iter_times, conv_1, conv_2 = f.readline().split()
            msg_len, repeat_times, iter_times, conv_1, conv_2 = int(msg_len), int(repeat_times), int(iter_times), int(conv_1), int(conv_2)
            # 从第二行开始是snr和ber
            for line in f.readlines():
                snr.append(float(line.split()[0]))
                ber.append(float(line.split()[1]))
        plt.semilogy(snr, ber, '-o', label='turbo_iter_'+str(iter_time))
    
    infostr = f'''
    信息序列长度: {msg_len}
    重复次数: {repeat_times}
    turbo译码迭代次数: {iter_times}
    卷积码类型: ({conv_1}, {conv_2})
    '''

    plt.text(0.6, 0.2, infostr, fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel('SNR(dB)')
    plt.ylabel('BER')
    plt.grid()
    plt.legend()
    plt.savefig('assets/figure_turbo.png')
