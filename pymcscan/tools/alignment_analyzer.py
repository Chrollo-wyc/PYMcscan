import os
from typing import List, Dict, Set
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class GeneFeatSimple:
    """简化的基因特征"""
    name: str
    mol: str
    mid: int
    in_blocks: int = 0
    cr_blocks: int = 0
    sp: Set[str] = None

    def __init__(self):
        self.sp = set()


class AlignmentAnalyzer:
    """多重比对分析器"""

    def __init__(self):
        self.gene_map: Dict[str, GeneFeatSimple] = {}
        self.allg: List[GeneFeatSimple] = []
        self.stat1: Dict[str, List[int]] = defaultdict(list)
        self.stat2: Dict[str, List[int]] = defaultdict(list)

    def read_gff(self, filepath: str):
        """读取GFF文件"""
        if not os.path.exists(filepath):
            print(f"Cannot read gff_file: {filepath}")
            exit(1)

        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    break

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                mol = parts[0]
                name = parts[1]
                try:
                    mid = int(parts[2])
                except ValueError:
                    continue

                spt = mol[:2]
                if spt not in self.stat1:
                    self.stat1[spt] = [0]
                    self.stat2[spt] = [0]

                gene = GeneFeatSimple()
                gene.name = name
                gene.mol = mol
                gene.mid = mid
                gene.in_blocks = 0
                gene.cr_blocks = 0

                self.gene_map[name] = gene

        # 创建排序的基因列表
        self.allg = sorted(self.gene_map.values(),
                           key=lambda x: (x.mol, x.mid))

    def add_block(self, s: GeneFeatSimple, t: GeneFeatSimple,
                  block_type: int, sp: str = ""):
        """添加块"""
        # 找到s和t之间的所有基因
        for gene in self.allg:
            if gene.mol == s.mol and s.mid <= gene.mid <= t.mid:
                if block_type == 0:
                    gene.in_blocks += 1
                else:
                    gene.cr_blocks += 1
                    if sp:
                        gene.sp.add(sp)

    def read_synteny(self, filepath: str):
        """读取共线性文件"""
        if not os.path.exists(filepath):
            print(f"Cannot read collinearity_file: {filepath}")
            exit(1)

        with open(filepath, 'r') as f:
            lines = f.readlines()

        i = 0
        strand = 1
        s1 = t1 = s2 = t2 = ""

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if line.startswith('## Alignment'):
                # 处理上一个块
                if i > 0 and s1 and t1 and s2 and t2:
                    gfs1 = self.gene_map.get(s1)
                    gft1 = self.gene_map.get(t1)
                    gfs2 = self.gene_map.get(s2)
                    gft2 = self.gene_map.get(t2)

                    if gfs1 and gft1 and gfs2 and gft2:
                        sp1 = gfs1.mol[:2]
                        sp2 = gfs2.mol[:2]
                        block_type = 0 if sp1 == sp2 else 1

                        self.add_block(gfs1, gft1, block_type, sp2)
                        self.add_block(gfs2, gft2, block_type, sp1)

                i += 1
                strand = 1
                if 'minus' in line:
                    strand = 0
                continue

            if line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                continue

            gene1 = parts[1]
            gene2 = parts[2]

            # 确定起止基因
            if s1 == "":
                s1 = gene1
                if strand == 1:
                    s2 = gene2
                else:
                    t2 = gene2
            else:
                t1 = gene1
                if strand == 1:
                    t2 = gene2
                else:
                    s2 = gene2

        # 处理最后一个块
        if s1 and t1 and s2 and t2:
            gfs1 = self.gene_map.get(s1)
            gft1 = self.gene_map.get(t1)
            gfs2 = self.gene_map.get(s2)
            gft2 = self.gene_map.get(t2)

            if gfs1 and gft1 and gfs2 and gft2:
                sp1 = gfs1.mol[:2]
                sp2 = gfs2.mol[:2]
                block_type = 0 if sp1 == sp2 else 1

                self.add_block(gfs1, gft1, block_type, sp2)
                self.add_block(gfs2, gft2, block_type, sp1)

    def print_file(self, filepath: str):
        """输出结果文件"""
        with open(filepath, 'w') as result:
            for gene in self.allg:
                # 更新统计
                spt = gene.mol[:2]

                # 更新self.stat1
                while len(self.stat1[spt]) <= gene.in_blocks:
                    self.stat1[spt].append(0)
                self.stat1[spt][gene.in_blocks] += 1

                # 更新self.stat2
                while len(self.stat2[spt]) <= gene.cr_blocks:
                    self.stat2[spt].append(0)
                self.stat2[spt][gene.cr_blocks] += 1

                # 输出基因信息
                sp_str = ",".join(gene.sp) if gene.sp else ""
                result.write(f"{gene.mol}\t{gene.name}\t{gene.in_blocks}\t{gene.cr_blocks}\t{len(gene.sp)}\t{sp_str}\n")

        # 输出统计信息
        print("Self-genome comparison:")
        print("Reference genome\tDuplication depth:gene number")
        for spt, counts in self.stat1.items():
            line = spt
            for depth, num in enumerate(counts):
                if num > 0:
                    line += f"\t{depth}:{num}"
            print(line)

        print("\nCross-genome comparison:")
        print("Reference genome\tDuplication depth:gene number")
        for spt, counts in self.stat2.items():
            line = spt
            for depth, num in enumerate(counts):
                if num > 0:
                    line += f"\t{depth}:{num}"
            print(line)

    def run(self, gff_file: str, synteny_file: str, output_file: str):
        """运行主流程"""
        self.read_gff(gff_file)
        self.read_synteny(synteny_file)
        self.print_file(output_file)