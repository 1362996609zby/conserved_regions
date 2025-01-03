from Bio import AlignIO
# conda activate hmm_env
# python conserved_regions.py
# 读取多序列比对文件
alignment_file = "BaiP_after_cdhit_aligned.fasta"
alignment = AlignIO.read(alignment_file, "fasta")

# 设置保守性阈值（如一致性比例90%）
threshold = 0.6

# 提取保守列
def get_conserved_columns(alignment, threshold):
    conserved_columns = []
    for col_index in range(alignment.get_alignment_length()):
        column = [seq[col_index] for seq in alignment]
        if column.count('-') / len(column) < (1 - threshold):  # 检查缺失情况
            most_common = max(set(column), key=column.count)
            if column.count(most_common) / len(column) >= threshold:
                conserved_columns.append(col_index)
    return conserved_columns

# 获取保守列
conserved_columns = get_conserved_columns(alignment, threshold)

# 提取保守序列
with open("BaiP_conserved_regions.fasta", "w") as output_file:
    for record in alignment:
        conserved_seq = "".join([record.seq[i] for i in conserved_columns])
        output_file.write(f">{record.id}\n{conserved_seq}\n")

print(f"保守区域已提取，共有{len(conserved_columns)}列")
