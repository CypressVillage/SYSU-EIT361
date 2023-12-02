import os

assets_path = 'assets'
build_path = 'build'

calc_actions = (
    'base',       # 计算并保存数据
    'turbo',      # 计算并保存不同迭代次数turbo译码数据
    'viterbi',    # 计算并保存不同结构卷积码编码数据
)
plot_actions = (
    'base',       # 画图
    'turbo',      # 画不同迭代次数turbo译码图
    'viterbi'     # 画不同结构卷积码编码图
)

arg_map = {
    'base': (
        'viterbi_hard',
        'viterbi_soft',
        'bcjr',
        'turbo'
    ),
    'turbo': (
        f'turbo {iter_time}'
        for iter_time in range(1, 10)
    ),
    'viterbi': (
        f'{method} {conv1} {conv2}'
        for method in (
            'viterbi_hard',
            'viterbi_soft'
        )
        for (conv1, conv2) in (
            (7, 5),
            (15, 13),
            (23, 35),
            (133, 171),
        )
    )
}
label_map = {
    'base': lambda arg: arg,
    'turbo': lambda arg: arg.replace('turbo', 'iter'),
    'viterbi': lambda arg: arg.replace('viterbi_', '').replace('hard', '硬判决').replace('soft', '软判决')
}


def path_with_os(path):
    is_linux = os.name == 'posix'
    if is_linux:
        return path.replace('.exe', '')
    else:
        return path.replace('/', '\\')

if not os.path.exists(assets_path):
    os.mkdir(assets_path)
if not os.path.exists(path_with_os(f'./{build_path}/main.exe')):
    print('请先编译项目')
    exit(0)


for action in calc_actions:
    for arg in arg_map[action]:
        output_file_name = path_with_os(f'{assets_path}/data_{action}_{arg}.txt'.replace(' ', '_'))
        # 删除旧文件
        if os.path.exists(output_file_name):
            os.remove(output_file_name)

        os.system(path_with_os(f'cd {assets_path} && ../{build_path}/main.exe {arg}'))
        os.rename(f'{assets_path}/data.txt', output_file_name)


if plot_actions:
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

for action in plot_actions:
    plt.figure()

    for arg in arg_map[action]:
        snr, ber = [], []
        output_file_name = path_with_os(f'{assets_path}/data_{action}_{arg}.txt'.replace(' ', '_'))
        with open(output_file_name) as f:
            # 第一行是message_length，重复次数，turbo译码迭代次数，卷积码类型
            msg_len, repeat_times, iter_times, conv_1, conv_2 = map(int, f.readline().split())
            # 从第二行开始是snr和ber
            for line in f.readlines():
                snr.append(float(line.split()[0]))
                ber.append(float(line.split()[1]))
        plt.semilogy(snr, ber, '-o', label=label_map[action](arg))

    infostr = f'''
    信息序列长度: {msg_len}
    重复次数: {repeat_times}
    turbo译码迭代次数: {iter_times}
    卷积码类型: ({conv_1}, {conv_2})
    '''
    plt.text(0.5, 0.5, infostr, fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel('SNR(dB)')
    plt.ylabel('BER')
    plt.grid()
    plt.legend()
    plt.savefig(f'{assets_path}/figure_{action}.png')
