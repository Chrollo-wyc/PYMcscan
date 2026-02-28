from pymcscan.core.data_structures import *


class GeneClassifier:
    """基因分类器"""

    def __init__(self, context: MCScanContext):
        self.context = context

    def cmpt_duptype(self):
        """计算重复类型"""
        # 初始化gene_more
        self.context.gene_more = [MoreFeature() for _ in range(len(self.context.allg))]

        # 第一步：所有匹配的基因标记为分散重复(1)
        for match in self.context.match_list:
            if match.gene1 in self.context.gene_map and match.gene2 in self.context.gene_map:
                gf1 = self.context.gene_map[match.gene1]
                gf2 = self.context.gene_map[match.gene2]

                idx1 = gf1.gene_id
                idx2 = gf2.gene_id

                if idx1 < len(self.context.gene_more) and self.context.gene_more[idx1].tandem < 1:
                    self.context.gene_more[idx1].tandem = 1
                if idx2 < len(self.context.gene_more) and self.context.gene_more[idx2].tandem < 1:
                    self.context.gene_more[idx2].tandem = 1

                # 如果是近端重复(2)
                if (abs(gf1.gene_id - gf2.gene_id) < self.context.N_PROXIMAL and
                        gf1.mol == gf2.mol):
                    if idx1 < len(self.context.gene_more) and self.context.gene_more[idx1].tandem < 2:
                        self.context.gene_more[idx1].tandem = 2
                    if idx2 < len(self.context.gene_more) and self.context.gene_more[idx2].tandem < 2:
                        self.context.gene_more[idx2].tandem = 2

                    # 如果是串联重复(3)
                    if abs(gf1.gene_id - gf2.gene_id) == 1:
                        if idx1 < len(self.context.gene_more):
                            self.context.gene_more[idx1].tandem = 3
                        if idx2 < len(self.context.gene_more):
                            self.context.gene_more[idx2].tandem = 3

        # 第二步：在共线性块中的基因标记为WGD/片段重复(4)
        for seg in self.context.seg_list:
            for pid in seg.pids:
                if pid < len(self.context.match_list):
                    match = self.context.match_list[pid]
                    if match.gene1 in self.context.gene_map:
                        gf1 = self.context.gene_map[match.gene1]
                        idx1 = gf1.gene_id
                        if idx1 < len(self.context.gene_more):
                            self.context.gene_more[idx1].tandem = 4

                    if match.gene2 in self.context.gene_map:
                        gf2 = self.context.gene_map[match.gene2]
                        idx2 = gf2.gene_id
                        if idx2 < len(self.context.gene_more):
                            self.context.gene_more[idx2].tandem = 4

    def print_cls(self, prefix_fn: str):
        """输出分类结果"""
        stat_fn = f"{prefix_fn}.gene_type"
        num = [0, 0, 0, 0, 0]  # 0-4类型

        with open(stat_fn, 'w') as result:
            for i, gene in enumerate(self.context.allg):
                tandem_type = 0
                if i < len(self.context.gene_more):
                    tandem_type = self.context.gene_more[i].tandem

                result.write(f"{gene.name}\t{tandem_type}\n")
                if tandem_type < len(num):
                    num[tandem_type] += 1

        print("Type of dup\tCode\tNumber")
        print(f"Singleton\t0\t{num[0]}")
        print(f"Dispersed\t1\t{num[1]}")
        print(f"Proximal\t2\t{num[2]}")
        print(f"Tandem\t3\t{num[3]}")
        print(f"WGD or segmental\t4\t{num[4]}")

    def cls_main(self, prefix_fn: str):
        """分类器主函数"""
        self.cmpt_duptype()
        self.print_cls(prefix_fn)