import os
from pymcscan.core.data_structures import *


class DataReader:
    """数据读取器"""

    def __init__(self, context: MCScanContext):
        self.context = context

    def read_gff(self, filepath: str):
        """读取GFF文件"""
        if not os.path.exists(filepath):
            print(f"Cannot read gff_file: {filepath}")
            exit(1)

        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                mol = parts[0]
                name = parts[1]
                try:
                    mid = int(parts[2])
                except ValueError:
                    continue

                gene = GeneFeature(name=name, mol=mol, mid=mid)
                self.context.gene_map[name] = gene

        # 填充allg并设置gene_id
        self.context.fill_allg()

    def read_blast(self, filepath: str):
        """读取BLAST文件"""
        if not os.path.exists(filepath):
            print(f"Cannot read blast_file: {filepath}")
            exit(1)

        blast_map = {}  # gene_pair -> best_evalue
        total_num = 0

        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 11:
                    continue

                gene1, gene2 = parts[0], parts[1]
                try:
                    evalue = float(parts[10])
                except (ValueError, IndexError):
                    continue

                # 跳过自比对
                if gene1 == gene2:
                    continue

                # 生成排序的基因对key
                if gene1 < gene2:
                    gene_key = f"{gene1}&{gene2}"
                else:
                    gene_key = f"{gene2}&{gene1}"

                # 保留最佳e值
                if gene_key not in blast_map or evalue < blast_map[gene_key]:
                    blast_map[gene_key] = evalue

                total_num += 1

        # 转换为BlastRecord
        pair_id = 0
        for gene_key, evalue in blast_map.items():
            gene1, gene2 = gene_key.split('&')

            # 检查基因是否在gene_map中
            if gene1 not in self.context.gene_map or gene2 not in self.context.gene_map:
                continue

            gf1 = self.context.gene_map[gene1]
            gf2 = self.context.gene_map[gene2]

            # 过滤条件
            if not gf1.mol or not gf2.mol:
                continue

            if self.context.IN_SYNTENY == 1 and gf1.mol[:2] != gf2.mol[:2]:
                continue

            # 确定基因顺序 - 修复字符串比较
            if gf1.mol < gf2.mol:
                # gf1.mol 在字典序上小于 gf2.mol
                mol_pair = f"{gf1.mol}&{gf2.mol}"
                rec = BlastRecord(gene1=gene1, gene2=gene2,
                                  mol_pair=mol_pair, pair_id=pair_id, score=evalue)
            elif gf1.mol == gf2.mol:
                # 相同染色体/分子
                if gf1.mid <= gf2.mid:
                    rec = BlastRecord(gene1=gene1, gene2=gene2,
                                      mol_pair=f"{gf1.mol}&{gf2.mol}",
                                      pair_id=pair_id, score=evalue)
                else:
                    rec = BlastRecord(gene1=gene2, gene2=gene1,
                                      mol_pair=f"{gf2.mol}&{gf1.mol}",
                                      pair_id=pair_id, score=evalue)
            else:
                # gf1.mol 在字典序上大于 gf2.mol
                rec = BlastRecord(gene1=gene2, gene2=gene1,
                                  mol_pair=f"{gf2.mol}&{gf1.mol}",
                                  pair_id=pair_id, score=evalue)

            # 更新mol_pairs计数
            if self.context.IN_SYNTENY != 2 or gf1.mol[:2] != gf2.mol[:2]:
                self.context.mol_pairs[rec.mol_pair] += 1

            self.context.match_list.append(rec)
            pair_id += 1

        selected_num = len(self.context.match_list)
        print(f"{selected_num} matches imported ({total_num - selected_num} discarded)")

    def read_orthomcl(self, filepath: str):
        """读取OrthoMCL同源文件"""
        # 实现类似read_blast，但支持不同的e_value模式
        pass