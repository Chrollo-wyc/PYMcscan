#!/usr/bin/env python3
"""
mycircos.draw - 基因组 Circos 图绘制模块
"""

import argparse
import os
import sys
import collections
import matplotlib.pyplot as plt
from pycircos import Garc, Gcircle



def parse_gff(gff_file):
    """
    解析 GFF 文件，提取染色体信息和基因位置

    GFF 格式（三列制表符分隔）:
        染色体ID    基因名称    位置

    返回:
        chromosomes: dict, 键为染色体ID，值为该染色体长度（最大位置）
        gene_positions: dict, 键为基因名称，值为 (染色体ID, 位置)
    """
    chromosomes = {}
    gene_positions = {}

    try:
        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    print(f"警告: 第 {line_num} 行格式错误，跳过: {line}", file=sys.stderr)
                    continue

                chrom = parts[0].strip()
                gene = parts[1].strip()
                try:
                    pos = int(parts[2])
                except ValueError:
                    print(f"警告: 第 {line_num} 行位置不是整数，跳过: {parts[2]}", file=sys.stderr)
                    continue

                # 记录基因位置
                gene_positions[gene] = (chrom, pos)

                # 更新染色体长度
                if chrom not in chromosomes or pos > chromosomes[chrom]:
                    chromosomes[chrom] = pos

    except Exception as e:
        print(f"解析 GFF 文件时出错: {e}", file=sys.stderr)
        raise

    return chromosomes, gene_positions


def parse_collinearity(coll_file, gene_positions):
    """
    解析 MCScanX 共线性文件，提取共线性块之间的连接信息

    共线性文件格式:
        ## Alignment 0: score=xxx ...
        i-j:    gene1    gene2    score

    返回:
        links: list, 每个元素为 (chrom1, start1, end1, chrom2, start2, end2)
    """
    links = []
    current_block_genes1 = []
    current_block_genes2 = []

    try:
        with open(coll_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    if line.startswith('## Alignment') and current_block_genes1:
                        # 处理上一个块
                        if len(current_block_genes1) >= 2 and len(current_block_genes2) >= 2:
                            # 使用块的首尾基因创建连接
                            g1_first, g1_last = current_block_genes1[0], current_block_genes1[-1]
                            g2_first, g2_last = current_block_genes2[0], current_block_genes2[-1]

                            if g1_first in gene_positions and g1_last in gene_positions and \
                                    g2_first in gene_positions and g2_last in gene_positions:
                                chrom1, pos1_first = gene_positions[g1_first]
                                _, pos1_last = gene_positions[g1_last]
                                chrom2, pos2_first = gene_positions[g2_first]
                                _, pos2_last = gene_positions[g2_last]

                                # 确保起始 < 结束
                                if pos1_first > pos1_last:
                                    pos1_first, pos1_last = pos1_last, pos1_first
                                if pos2_first > pos2_last:
                                    pos2_first, pos2_last = pos2_last, pos2_first

                                links.append((chrom1, pos1_first, pos1_last,
                                              chrom2, pos2_first, pos2_last))

                        current_block_genes1 = []
                        current_block_genes2 = []
                    continue

                # 解析比对行
                if '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        # 格式: i-j:    gene1    gene2    score
                        gene_part = parts[1].strip()
                        if ' ' in gene_part:
                            gene1, gene2 = gene_part.split()[:2]
                        else:
                            # 备用格式
                            continue

                        if gene1 in gene_positions and gene2 in gene_positions:
                            chrom1, _ = gene_positions[gene1]
                            chrom2, _ = gene_positions[gene2]

                            # 按染色体分组存储
                            if not current_block_genes1 or chrom1 == gene_positions[current_block_genes1[-1]][0]:
                                current_block_genes1.append(gene1)
                                current_block_genes2.append(gene2)
                            else:
                                # 新块开始
                                if len(current_block_genes1) >= 2:
                                    g1_first, g1_last = current_block_genes1[0], current_block_genes1[-1]
                                    g2_first, g2_last = current_block_genes2[0], current_block_genes2[-1]

                                    if g1_first in gene_positions and g1_last in gene_positions:
                                        chrom1, pos1_first = gene_positions[g1_first]
                                        _, pos1_last = gene_positions[g1_last]
                                        chrom2, pos2_first = gene_positions[g2_first]
                                        _, pos2_last = gene_positions[g2_last]

                                        if pos1_first > pos1_last:
                                            pos1_first, pos1_last = pos1_last, pos1_first
                                        if pos2_first > pos2_last:
                                            pos2_first, pos2_last = pos2_last, pos2_first

                                        links.append((chrom1, pos1_first, pos1_last,
                                                      chrom2, pos2_first, pos2_last))

                                current_block_genes1 = [gene1]
                                current_block_genes2 = [gene2]

        # 处理最后一个块
        if len(current_block_genes1) >= 2 and len(current_block_genes2) >= 2:
            g1_first, g1_last = current_block_genes1[0], current_block_genes1[-1]
            g2_first, g2_last = current_block_genes2[0], current_block_genes2[-1]

            if g1_first in gene_positions and g1_last in gene_positions:
                chrom1, pos1_first = gene_positions[g1_first]
                _, pos1_last = gene_positions[g1_last]
                chrom2, pos2_first = gene_positions[g2_first]
                _, pos2_last = gene_positions[g2_last]

                if pos1_first > pos1_last:
                    pos1_first, pos1_last = pos1_last, pos1_first
                if pos2_first > pos2_last:
                    pos2_first, pos2_last = pos2_last, pos2_first

                links.append((chrom1, pos1_first, pos1_last,
                              chrom2, pos2_first, pos2_last))

    except Exception as e:
        print(f"解析共线性文件时出错: {e}", file=sys.stderr)
        raise

    return links


# ---------- 核心绘图函数（新添加）----------
def draw_circos(gff_file, collinearity_file, output_file, title="Genome Circos Plot",
                figsize=(10,10), dpi=300, link_color="red", link_alpha=0.3, show_labels=True):
    """
    使用 mycircos 绘制 Circos 图
    :param gff_file: GFF文件路径
    :param collinearity_file: 共线性文件路径
    :param output_file: 输出图片路径
    :param title: 图表标题
    :param figsize: 图片尺寸 (宽, 高) 英寸
    :param dpi: 图片分辨率
    :param link_color: 连接线颜色
    :param link_alpha: 连接线透明度
    :param show_labels: 是否显示染色体标签
    """
    # 解析文件
    chromosomes, gene_positions = parse_gff(gff_file)
    links = parse_collinearity(collinearity_file, gene_positions)

    # 创建画布
    circle = Gcircle(figsize=figsize)

    # 添加染色体弧
    for chrom_name, chrom_size in chromosomes.items():
        arc = Garc(
            arc_id=chrom_name,
            size=chrom_size,
            interspace=2,
            raxis_range=(800, 950),
            labelposition=70,
            label_visible=show_labels
        )
        circle.add_garc(arc)

    # 设置角度
    circle.set_garcs(-65, 245)

    # 添加刻度轨道
    for chrom_name in chromosomes:
        circle.tickplot(
            chrom_name,
            raxis_range=(950, 980),
            tickinterval=chromosomes[chrom_name] // 5,
            ticklabels=None
        )

    # 绘制连接线
    for link in links:
        chrom1, start1, end1, chrom2, start2, end2 = link
        try:
            source = (chrom1, start1, end1, 615)
            destination = (chrom2, start2, end2, 615)
            circle.chord_plot(source, destination,
                              facecolor=link_color, alpha=link_alpha)
        except Exception as e:
            print(f"警告: 绘制连接线失败 {e}", file=sys.stderr)

    # 设置标题
    plt.title(title, fontsize=14, y=1.08)

    # 保存
    circle.figure.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"Circos 图已保存至: {output_file}")

# ---------- 命令行入口（保持不变）----------
def main():
    parser = argparse.ArgumentParser(description="绘制基因组 Circos 图")
    parser.add_argument("--gff", required=True, help="GFF 文件路径")
    parser.add_argument("--collinearity", required=True, help="共线性文件路径")
    parser.add_argument("--output", required=True, help="输出图片路径")
    parser.add_argument("--title", default="Genome Circos Plot", help="图表标题")
    parser.add_argument("--figwidth", type=float, default=10, help="图片宽度 (英寸)")
    parser.add_argument("--figheight", type=float, default=10, help="图片高度 (英寸)")
    parser.add_argument("--dpi", type=int, default=300, help="图片 DPI")
    parser.add_argument("--link-color", default="red", help="连接线颜色")
    parser.add_argument("--link-alpha", type=float, default=0.3, help="连接线透明度")
    args = parser.parse_args()

    draw_circos(
        gff_file=args.gff,
        collinearity_file=args.collinearity,
        output_file=args.output,
        title=args.title,
        figsize=(args.figwidth, args.figheight),
        dpi=args.dpi,
        link_color=args.link_color,
        link_alpha=args.link_alpha,
        show_labels=True  # 默认显示标签
    )

if __name__ == "__main__":
    main()