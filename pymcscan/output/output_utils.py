# [file name]: mcscanx/output/output_utils.py
"""
输出工具模块
包含比对结果输出、参数输出等功能。
"""

import os
from typing import Set, TextIO
from pymcscan.core.data_structures import MCScanContext, SegFeature
from pymcscan.utils.file_utils import progress


class OutputUtils:
    """输出工具类"""

    def __init__(self, context: MCScanContext):
        self.context = context

    def print_params(self, f: TextIO) -> None:
        """输出算法参数"""
        f.write("############### Parameters ###############\n")
        f.write(f"# MATCH_SCORE: {self.context.MATCH_SCORE}\n")
        f.write(f"# MATCH_SIZE: {self.context.MATCH_SIZE}\n")
        f.write(f"# GAP_PENALTY: {self.context.GAP_PENALTY}\n")
        f.write(f"# OVERLAP_WINDOW: {self.context.OVERLAP_WINDOW}\n")
        f.write(f"# E_VALUE: {self.context.E_VALUE}\n")
        f.write(f"# MAX GAPS: {self.context.MAX_GAPS}\n")
        f.write("##########################################\n\n")

    def print_align(self, output_file: str) -> None:
        """输出比对结果到文件"""
        with open(output_file, 'w') as f:
            self.print_params(f)

            # 统计共线性基因
            colgenes: Set[str] = set()
            for seg in self.context.seg_list:
                for pid in seg.pids:
                    if pid < len(self.context.match_list):
                        match = self.context.match_list[pid]
                        colgenes.add(match.gene1)
                        colgenes.add(match.gene2)

            f.write("############### Statistics ###############\n")
            if self.context.gene_map:
                total_genes = len(self.context.gene_map)
                colinear_count = len(colgenes)
                percentage = 100.0 * colinear_count / total_genes if total_genes > 0 else 0.0
                f.write(f"# Number of collinear genes: {colinear_count}, Percentage: {percentage:.2f}\n")
                f.write(f"# Number of all genes: {total_genes}\n")

            # 如果是OrthoMCL版本，输出物种间统计
            if hasattr(self.context, 'cmp_sp') and self.context.cmp_sp:
                f.write("\n# Cross-species comparison statistics:\n")
                for sp_pair, stat in self.context.cmp_sp.items():
                    if stat.all_num > 0:
                        percentage = 100.0 * stat.syn_num / stat.all_num
                        f.write(f"# {sp_pair}: {stat.syn_num} / {stat.all_num} ({percentage:.2f}%)\n")

            f.write("##########################################\n\n")

            # 输出每个共线性块
            for i, seg in enumerate(self.context.seg_list):
                nanchor = len(seg.pids)
                strand_str = "plus" if seg.same_strand else "minus"

                f.write(f"## Alignment {i}: score={seg.score:.1f} e_value={seg.e_value:.2g} ")
                f.write(f"N={nanchor} {seg.mol_pair} {strand_str}\n")

                for j, pid in enumerate(seg.pids):
                    if pid < len(self.context.match_list):
                        match = self.context.match_list[pid]
                        f.write(f"{i:3d}-{j:3d}:\t{match.gene1}\t{match.gene2}\t{match.score:7.1g}\n")

    def print_gene_types(self, output_file: str) -> None:
        """输出基因类型结果"""
        if not self.context.gene_more or len(self.context.allg) != len(self.context.gene_more):
            err("Gene classification data not properly initialized")
            return

        num = [0, 0, 0, 0, 0]  # 0-4类型

        with open(output_file, 'w') as f:
            for i, gene in enumerate(self.context.allg):
                tandem_type = 0
                if i < len(self.context.gene_more):
                    tandem_type = self.context.gene_more[i].tandem

                f.write(f"{gene.name}\t{tandem_type}\n")
                if tandem_type < len(num):
                    num[tandem_type] += 1

        # 打印统计信息到控制台
        progress("Type of dup\tCode\tNumber")
        progress(f"Singleton\t0\t{num[0]}")
        progress(f"Dispersed\t1\t{num[1]}")
        progress(f"Proximal\t2\t{num[2]}")
        progress(f"Tandem\t3\t{num[3]}")
        progress(f"WGD or segmental\t4\t{num[4]}")

    def print_tandem_arrays(self, anchor1: str, array1: str, anchor2: str, array2: str,
                            output_file: str) -> None:
        """输出串联阵列结果"""
        with open(output_file, 'w') as f:
            f.write("Anchor1\tTandem_array1\tAnchor2\tTandem_array2\n")
            for a1, arr1, a2, arr2 in zip(anchor1, array1, anchor2, array2):
                f.write(f"{a1}\t{arr1}\t{a2}\t{arr2}\n")

    def print_html_header(self, f: TextIO, max_level: int, title: str = "MCScanX Results") -> None:
        """输出HTML文件头部"""
        f.write(f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>{title}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .tandem {{ color: #ff0000; }}
        .collinear {{ background-color: #ffff99; }}
        .legend {{ margin: 20px 0; padding: 10px; border: 1px solid #ccc; }}
    </style>
</head>
<body>
    <h1>{title}</h1>
    <div class="legend">
        <h3>Legend:</h3>
        <p><span style="color: #ee0000;">Red</span>: Tandem duplicates</p>
        <p><span style="background-color: #ffff99;">Yellow</span>: Collinear genes</p>
    </div>
    <table cellspacing='0' cellpadding='0' align='left'>
        <tr align='center'>
            <td>Duplication depth</td>
            <td>&nbsp;&nbsp;Reference chromosome</td>
            <td align='left' colspan='{2 * max_level}'>&nbsp;&nbsp;Collinear blocks</td>
        </tr>
""")

    def print_html_footer(self, f: TextIO) -> None:
        """输出HTML文件尾部"""
        f.write("""
    </table>
</body>
</html>
""")

    def generate_report(self, output_file: str, include_stats: bool = True) -> None:
        """生成综合报告"""
        with open(output_file, 'w') as f:
            f.write(f"MCScanX Analysis Report\n")
            f.write(f"Generated: {os.path.basename(output_file)}\n")
            f.write("=" * 60 + "\n\n")

            # 参数部分
            f.write("PARAMETERS:\n")
            f.write("-" * 40 + "\n")
            self.print_params(f)

            if include_stats:
                # 统计部分
                f.write("\nSTATISTICS:\n")
                f.write("-" * 40 + "\n")

                # 基因统计
                total_genes = len(self.context.gene_map)
                f.write(f"Total genes: {total_genes}\n")

                # 比对统计
                total_alignments = len(self.context.seg_list)
                f.write(f"Total alignments: {total_alignments}\n")

                # 分子对统计
                f.write(f"Molecular pairs analyzed: {len(self.context.mol_pairs)}\n")

                # 如果是分类器模式，输出分类统计
                if self.context.gene_more and len(self.context.allg) == len(self.context.gene_more):
                    f.write("\nDUPLICATION CLASSIFICATION:\n")
                    types = ["Singleton", "Dispersed", "Proximal", "Tandem", "WGD/Segmental"]
                    counts = [0] * 5

                    for i, gene in enumerate(self.context.allg):
                        tandem_type = self.context.gene_more[i].tandem if i < len(self.context.gene_more) else 0
                        if 0 <= tandem_type < 5:
                            counts[tandem_type] += 1

                    for i, (type_name, count) in enumerate(zip(types, counts)):
                        f.write(f"{type_name}: {count} ({100.0 * count / total_genes:.2f}%)\n")