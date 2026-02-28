#!/usr/bin/env python3
"""
重复基因分类器命令行工具
用于分类重复基因。
用法: python -m pymcscan.cli.dup_classifier prefix_fn [options]
"""

import argparse
import sys
import os

# 添加项目根目录到路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from pymcscan.core.data_structures import MCScanContext, ScoreTuple
from pymcscan.core.reader import DataReader
from pymcscan.core.dagchainer import DAGChainer
from pymcscan.core.classifier import GeneClassifier
from pymcscan.utils.file_utils import progress, uglyTime, errAbort


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="MCScanX Duplicate Gene Classifier",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
参数说明:
  -k MATCH_SCORE     最终分数=MATCH_SCORE+NUM_GAPS*GAP_PENALTY (默认: 50)
  -g GAP_PENALTY     缺口惩罚 (默认: -1)
  -s MATCH_SIZE      调用共线性块所需的基因数 (默认: 5)
  -e E_VALUE         比对显著性 (默认: 1e-5)
  -w OVERLAP_WINDOW  折叠BLAST匹配的最大距离 (默认: 5)
  -n N_PROXIMAL      调用近端重复的最大距离 (默认: 10)
  -m MAX_GAPS        允许的最大缺口数 (默认: 25)
        """
    )

    parser.add_argument("prefix_fn", help="输入文件前缀（不含扩展名）")
    parser.add_argument("-k", "--match-score", type=int, default=50,
                        help="匹配分数")
    parser.add_argument("-g", "--gap-penalty", type=int, default=-1,
                        help="缺口惩罚")
    parser.add_argument("-s", "--match-size", type=int, default=5,
                        help="匹配大小")
    parser.add_argument("-e", "--evalue", type=float, default=1e-5,
                        help="E值")
    parser.add_argument("-w", "--overlap-window", type=int, default=5,
                        help="重叠窗口")
    parser.add_argument("-n", "--n-proximal", type=int, default=10,
                        help="近端距离")
    parser.add_argument("-m", "--max-gaps", type=int, default=25,
                        help="最大缺口数")

    return parser.parse_args()


def main():
    """主函数"""
    # 开始计时
    uglyTime(None)

    args = parse_arguments()

    # 初始化上下文
    context = MCScanContext()

    # 设置参数
    context.MATCH_SCORE = args.match_score
    context.MATCH_SIZE = args.match_size
    context.GAP_PENALTY = args.gap_penalty
    context.E_VALUE = args.evalue
    context.OVERLAP_WINDOW = args.overlap_window
    context.MAX_GAPS = args.max_gaps
    context.N_PROXIMAL = args.n_proximal
    context.set_cutoff()

    # 创建组件
    reader = DataReader(context)
    dagchainer = DAGChainer(context)
    classifier = GeneClassifier(context)

    # 读取数据
    progress("Reading GFF file...")
    gff_file = f"{args.prefix_fn}.gff"
    if not os.path.exists(gff_file):
        errAbort("GFF file not found: %s", gff_file)
    reader.read_gff(gff_file)

    progress("Reading BLAST file...")
    blast_file = f"{args.prefix_fn}.blast"
    if not os.path.exists(blast_file):
        errAbort("BLAST file not found: %s", blast_file)
    reader.read_blast(blast_file)

    progress("%d pairwise comparisons", len(context.mol_pairs))

    # 运行DAG链式算法
    for mol_pair, count in context.mol_pairs.items():
        if count >= context.MATCH_SIZE:
            # 构建score列表
            score_list = []
            for match in context.match_list:
                if match.mol_pair == mol_pair:
                    if match.gene1 in context.gene_map and match.gene2 in context.gene_map:
                        gf1 = context.gene_map[match.gene1]
                        gf2 = context.gene_map[match.gene2]

                        score_tuple = ScoreTuple(
                            pair_id=match.pair_id,
                            x=gf1.gene_id,
                            y=gf2.gene_id,
                            x_mid=gf1.mid,
                            y_mid=gf2.mid,
                            score=context.MATCH_SCORE,
                            gene1=match.gene1,
                            gene2=match.gene2
                        )
                        score_list.append(score_tuple)

            if score_list:
                dagchainer.dag_main(score_list, mol_pair)

    progress("%d alignments generated", len(context.seg_list))

    # 分类重复基因
    progress("Classifying duplicate genes...")
    classifier.cls_main(args.prefix_fn)

    uglyTime("Done!")


if __name__ == "__main__":
    main()