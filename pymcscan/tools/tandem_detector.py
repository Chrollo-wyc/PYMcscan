import os
from typing import List, Dict
from dataclasses import dataclass



@dataclass
class BlastRecordSimple:
    """简化的BLAST记录"""
    gene1: str
    gene2: str
    mol1: str
    mol2: str
    id1: int
    id2: int
    score: float


@dataclass
class TandemPair:
    """串联对"""
    gene1: str
    gene2: str
    sp: str
    mol: str
    id: int


@dataclass
class TandemArray:
    """串联阵列"""
    genes: List[str]
    mol: str
    sp: str


@dataclass
class TandemCluster:
    """串联簇"""
    anchor1: str
    anchor2: str
    array_id1: int = -1
    array_id2: int = -1


class TandemDetector:
    """串联阵列检测器"""

    def __init__(self):
        self.gene_map: Dict[str, Dict] = {}
        self.match_list: List[BlastRecordSimple] = []
        self.allg: List[Dict] = []

        # 串联相关数据结构
        self.alltandempair: List[TandemPair] = []
        self.alltandemarray: List[TandemArray] = []
        self.alltandemcluster: List[TandemCluster] = []
        self.tandengeneid: Dict[str, int] = {}

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

                self.gene_map[name] = {
                    'name': name,
                    'mol': mol,
                    'mid': mid,
                    'geneid': 0
                }

        # 创建allg并设置geneid
        self.allg = sorted(self.gene_map.values(),
                           key=lambda x: (x['mol'], x['mid']))

        for i, gene in enumerate(self.allg):
            gene['geneid'] = i
            self.gene_map[gene['name']]['geneid'] = i

    def read_blast(self, filepath: str):
        """读取BLAST文件"""
        if not os.path.exists(filepath):
            print(f"Cannot read blast_file: {filepath}")
            exit(1)

        blast_map = {}
        total_num = 0

        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    break

                parts = line.split('\t')
                if len(parts) < 11:
                    continue

                gene1, gene2 = parts[0], parts[1]
                try:
                    evalue = float(parts[10])
                except (ValueError, IndexError):
                    continue

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

        # 转换为记录
        for gene_key, evalue in blast_map.items():
            gene1, gene2 = gene_key.split('&')

            if gene1 not in self.gene_map or gene2 not in self.gene_map:
                continue

            gf1 = self.gene_map[gene1]
            gf2 = self.gene_map[gene2]

            if not gf1['mol'] or not gf2['mol']:
                continue

            # 确定顺序
            if gf1['geneid'] < gf2['geneid']:
                rec = BlastRecordSimple(
                    gene1=gene1, gene2=gene2,
                    mol1=gf1['mol'], mol2=gf2['mol'],
                    id1=gf1['geneid'], id2=gf2['geneid'],
                    score=evalue
                )
            else:
                rec = BlastRecordSimple(
                    gene1=gene2, gene2=gene1,
                    mol1=gf2['mol'], mol2=gf1['mol'],
                    id1=gf2['geneid'], id2=gf1['geneid'],
                    score=evalue
                )

            self.match_list.append(rec)

        selected_num = len(self.match_list)
        print(f"{selected_num} matches imported ({total_num - selected_num} discarded)")

    def compute_tandem(self):
        """计算串联阵列"""
        print("Detecting tandem arrays...")

        # 收集所有串联对
        for match in self.match_list:
            if (match.mol1 == match.mol2 and
                    (match.id1 - match.id2) == -1):
                pair = TandemPair(
                    gene1=match.gene1,
                    gene2=match.gene2,
                    mol=match.mol1,
                    sp=match.mol1[:2],
                    id=match.id1
                )
                self.alltandempair.append(pair)

        # 按mol和id排序
        self.alltandempair.sort(key=lambda x: (x.mol, x.id))

        # 构建串联阵列
        if not self.alltandempair:
            return

        temp_array = TandemArray(
            genes=[self.alltandempair[0].gene1, self.alltandempair[0].gene2],
            mol=self.alltandempair[0].mol,
            sp=self.alltandempair[0].sp
        )

        prev_id = self.alltandempair[0].id

        for i in range(1, len(self.alltandempair)):
            if self.alltandempair[i].id == prev_id + 1:
                temp_array.genes.append(self.alltandempair[i].gene2)
                prev_id += 1
            else:
                self.alltandemarray.append(temp_array)
                temp_array = TandemArray(
                    genes=[self.alltandempair[i].gene1, self.alltandempair[i].gene2],
                    mol=self.alltandempair[i].mol,
                    sp=self.alltandempair[i].sp
                )
                prev_id = self.alltandempair[i].id

        self.alltandemarray.append(temp_array)

        # 构建基因到阵列的映射
        for i, array in enumerate(self.alltandemarray):
            for gene in array.genes:
                self.tandengeneid[gene] = i

    def read_synteny(self, filepath: str):
        """读取共线性文件"""
        if not os.path.exists(filepath):
            print(f"Cannot read collinearity_file: {filepath}")
            exit(1)

        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                gene1 = parts[1]
                gene2 = parts[2]

                array_id1 = self.tandengeneid.get(gene1, -1)
                array_id2 = self.tandengeneid.get(gene2, -1)

                if array_id1 != -1 or array_id2 != -1:
                    cluster = TandemCluster(
                        anchor1=gene1,
                        anchor2=gene2,
                        array_id1=array_id1,
                        array_id2=array_id2
                    )
                    self.alltandemcluster.append(cluster)

    def print_file(self, filepath: str):
        """输出结果文件"""
        print("Outputting results...")

        with open(filepath, 'w') as result:
            result.write("Anchor1\tTandem_array1\tAnchor2\tTandem_array2\n")

            for cluster in self.alltandemcluster:
                result.write(f"{cluster.anchor1}\t")

                if cluster.array_id1 < 0:
                    result.write(f"{cluster.anchor1}\t")
                else:
                    genes = self.alltandemarray[cluster.array_id1].genes
                    result.write(",".join(genes) + "\t")

                result.write(f"{cluster.anchor2}\t")

                if cluster.array_id2 < 0:
                    result.write(f"{cluster.anchor2}\n")
                else:
                    genes = self.alltandemarray[cluster.array_id2].genes
                    result.write(",".join(genes) + "\n")

        print("Done!")

    def run(self, gff_file: str, blast_file: str, synteny_file: str, output_file: str):
        """运行主流程"""
        self.read_gff(gff_file)
        self.read_blast(blast_file)
        self.compute_tandem()
        self.read_synteny(synteny_file)
        self.print_file(output_file)