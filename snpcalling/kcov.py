# -*- coding: utf-8 -*-
import sys

# 检查命令行参数
if len(sys.argv) != 3:
    print(" python script.py input output")
    sys.exit(1)

# 获取输入和输出文件路径
input_file = sys.argv[1]
output_file = sys.argv[2]

# 读取文件前1000行并存储数据
data = []
with open(input_file, 'r') as file:
    for line_num, line in enumerate(file):
        if line_num >= 1000:
            break  # 只读取前1000行
        columns = line.split()
        data.append((int(columns[0]), int(columns[1])))

# 计算总和
total_sum = sum([point[1] for point in data])

# 初始化存储凸起部分的列表
peaks = []
cumulative_sum = data[0][1]

# 遍历数据，查找符合条件的凸起部分
# 每次检查连续三个点的趋势：先上升后下降
for i in range(1, len(data) - 1):
    # 更新累积和并检查是否超过95%
    cumulative_sum += data[i][1]
    if cumulative_sum / total_sum > 0.95:
        break  # 超过95%，停止查找，完全跳出循环，不再处理数据

    if data[i-1][1] < data[i][1] > data[i+1][1]:
        peaks.append(data[i][0])  # 记录凸起部分的横坐标
    if len(peaks) >= 2:  # 如果已经找到两个凸起，则停止查找
        break

# 检查第二个凸起是否满足条件：第二个凸起的横坐标/第一个凸起的横坐标在1.5-2.5之间
if len(peaks) == 2:
    ratio = peaks[1] / peaks[0]
    if ratio < 1.8 or ratio > 2.2:
        peaks.pop(1)  # 不符合条件，删去第二个凸起

# 将结果写入输出文件
with open(output_file, 'w') as outfile:
    if peaks:
        outfile.write(f"n:{len(peaks)}\n")
        for peak in peaks:
            outfile.write(f"{peak}\n")
    else:
        outfile.write("n:0\n")
        outfile.write("No peaks found! coverage too low!\n")
