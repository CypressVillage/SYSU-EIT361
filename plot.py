import matplotlib.pyplot as plt

# 读取data.txt文件中的SNR和BER数据
snr = []
ber = []
with open('data.txt', 'r') as f:
    for line in f.readlines():
        snr.append(float(line.split()[0]))
        ber.append(float(line.split()[1]))

# 绘制SNR-BER曲线
plt.figure(figsize=(8, 6))
plt.semilogy(snr, ber, 'b-o')
plt.xlabel('SNR(dB)')
plt.ylabel('BER')
plt.grid(True)
plt.show()
