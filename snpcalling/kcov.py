# -*- coding: utf-8 -*-
import sys

# ��������в���
if len(sys.argv) != 3:
    print(" python script.py input output")
    sys.exit(1)

# ��ȡ���������ļ�·��
input_file = sys.argv[1]
output_file = sys.argv[2]

# ��ȡ�ļ�ǰ1000�в��洢����
data = []
with open(input_file, 'r') as file:
    for line_num, line in enumerate(file):
        if line_num >= 1000:
            break  # ֻ��ȡǰ1000��
        columns = line.split()
        data.append((int(columns[0]), int(columns[1])))

# �����ܺ�
total_sum = sum([point[1] for point in data])

# ��ʼ���洢͹�𲿷ֵ��б�
peaks = []
cumulative_sum = data[0][1]

# �������ݣ����ҷ���������͹�𲿷�
# ÿ�μ����������������ƣ����������½�
for i in range(1, len(data) - 1):
    # �����ۻ��Ͳ�����Ƿ񳬹�95%
    cumulative_sum += data[i][1]
    if cumulative_sum / total_sum > 0.95:
        break  # ����95%��ֹͣ���ң���ȫ����ѭ�������ٴ�������

    if data[i-1][1] < data[i][1] > data[i+1][1]:
        peaks.append(data[i][0])  # ��¼͹�𲿷ֵĺ�����
    if len(peaks) >= 2:  # ����Ѿ��ҵ�����͹����ֹͣ����
        break

# ���ڶ���͹���Ƿ������������ڶ���͹��ĺ�����/��һ��͹��ĺ�������1.5-2.5֮��
if len(peaks) == 2:
    ratio = peaks[1] / peaks[0]
    if ratio < 1.8 or ratio > 2.2:
        peaks.pop(1)  # ������������ɾȥ�ڶ���͹��

# �����д������ļ�
with open(output_file, 'w') as outfile:
    if peaks:
        outfile.write(f"n:{len(peaks)}\n")
        for peak in peaks:
            outfile.write(f"{peak}\n")
    else:
        outfile.write("n:0\n")
        outfile.write("No peaks found! coverage too low!\n")
