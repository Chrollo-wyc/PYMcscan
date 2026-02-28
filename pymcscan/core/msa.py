import os
from pymcscan.core.data_structures import *


@dataclass
class NewEndpoint:
    """新端点"""
    n: GeneFeature
    index: int
    start: bool
    e: GeneFeature

    def __lt__(self, other: 'NewEndpoint') -> bool:
        return (self.n.mol < other.n.mol) or \
            (self.n.mol == other.n.mol and self.n.mid < other.n.mid)


class MSA:
    """多重序列比对"""

    def __init__(self, context: MCScanContext):
        self.context = context
        self.max_level = 1

    def get_endpoints(self):
        """获取端点"""
        endpoints = []
        n = len(self.context.seg_list)

        for i, seg in enumerate(self.context.seg_list):
            seg.index = i

            # 添加四个端点
            endpoints.append(NewEndpoint(
                n=seg.s1, index=2 * i, start=True, e=seg.t1
            ))
            endpoints.append(NewEndpoint(
                n=seg.t1, index=2 * i, start=False, e=seg.s1
            ))
            endpoints.append(NewEndpoint(
                n=seg.s2, index=2 * i + 1, start=True, e=seg.t2
            ))
            endpoints.append(NewEndpoint(
                n=seg.t2, index=2 * i + 1, start=False, e=seg.s2
            ))

        endpoints.sort()
        return endpoints

    def add_block(self, s: GeneFeature, t: GeneFeature, level: int):
        """添加块"""
        # 找到s和t之间的所有基因
        for gene in self.context.allg:
            if (gene.mol == s.mol and
                    s.mid <= gene.mid <= t.mid):
                # 确保cursor列表足够长
                while len(gene.cursor) < level:
                    gene.cursor.append(0)
                gene.cursor[level - 1] = 1

    def add_matchpoints(self, seg_index: int, level: int):
        """添加匹配点"""
        seg_idx = seg_index // 2
        if seg_idx >= len(self.context.seg_list):
            return

        seg = self.context.seg_list[seg_idx]

        if seg_index % 2 == 0:
            # 第一个基因组
            for pid in seg.pids:
                gene_name = self.context.match_list[pid].gene1
                if gene_name in self.context.gene_map:
                    gene = self.context.gene_map[gene_name]
                    if len(gene.cursor) >= level:
                        gene.cursor[level - 1] = pid + 2
        else:
            # 第二个基因组
            for pid in seg.pids:
                gene_name = self.context.match_list[pid].gene2
                if gene_name in self.context.gene_map:
                    gene = self.context.gene_map[gene_name]
                    if len(gene.cursor) >= level:
                        gene.cursor[level - 1] = -(pid + 2)

    def traverse(self, endpoints: List[NewEndpoint]):
        """遍历端点"""
        for endpoint in endpoints:
            if endpoint.start:
                gene = endpoint.n
                k = len(gene.cursor)

                if k == 0:
                    self.add_block(gene, endpoint.e, 1)
                    self.add_matchpoints(endpoint.index, 1)
                else:
                    # 找到第一个为0的位置
                    lev = next((j + 1 for j, val in enumerate(gene.cursor) if val == 0), k + 1)
                    self.add_block(gene, endpoint.e, lev)
                    self.add_matchpoints(endpoint.index, lev)

                    if lev > self.max_level:
                        self.max_level = lev

    def mark_tandem(self, prefix_fn: str):
        """标记串联重复"""
        # 初始化gene_more
        self.context.gene_more = [MoreFeature() for _ in range(len(self.context.allg))]

        # 计算深度
        for i, gene in enumerate(self.context.allg):
            depth = sum(1 for val in gene.cursor if val != 0)
            self.context.gene_more[i].depth = depth

        # 检测串联重复
        tandem_pairs = []
        for i, match in enumerate(self.context.match_list):
            if match.gene1 in self.context.gene_map and match.gene2 in self.context.gene_map:
                gf1 = self.context.gene_map[match.gene1]
                gf2 = self.context.gene_map[match.gene2]

                if (abs(gf1.gene_id - gf2.gene_id) == 1 and
                        gf1.mol == gf2.mol):
                    idx1 = gf1.gene_id
                    idx2 = gf2.gene_id

                    if idx1 < len(self.context.gene_more):
                        self.context.gene_more[idx1].tandem = 1
                    if idx2 < len(self.context.gene_more):
                        self.context.gene_more[idx2].tandem = 1

                    tandem_pairs.append((match.gene1, match.gene2))

        # 输出串联对
        if tandem_pairs:
            tandem_file = f"{prefix_fn}.tandem"
            print(f"Tandem pairs written to {tandem_file}")

            with open(tandem_file, 'w') as f:
                for gene1, gene2 in tandem_pairs:
                    f.write(f"{gene1},{gene2}\n")

    def print_html(self, prefix_fn: str):
        """输出HTML文件"""
        html_dir = f"{prefix_fn}.html"
        os.makedirs(html_dir, exist_ok=True)

        prev_mol = ""
        html_file = None

        for gene in self.context.allg:
            if gene.mol != prev_mol:
                if html_file:
                    html_file.write("</table></html>")
                    html_file.close()

                html_path = os.path.join(html_dir, f"{gene.mol}.html")
                print(f"Creating HTML file: {html_path}")

                html_file = open(html_path, 'w')
                html_file.write("<html><table cellspacing='0' cellpadding='0' align='left'>")
                html_file.write(
                    f"<tr align='center'><td>Duplication depth</td><td>&nbsp;&nbsp;Reference chromosome</td><td align='left' colspan='{2 * self.max_level}'>&nbsp;&nbsp;Collinear blocks</td></tr>\n")
                prev_mol = gene.mol

            # 确定颜色
            gene_idx = gene.gene_id
            color = "'#dddddd'"
            if gene_idx < len(self.context.gene_more) and self.context.gene_more[gene_idx].tandem:
                color = "'#ee0000'"

            # 开始行
            depth = self.context.gene_more[gene_idx].depth if gene_idx < len(self.context.gene_more) else 0
            html_file.write(f"<tr align='center'><td>{depth}</td><td bgcolor={color}>{gene.name}</td>")

            # 输出cursor信息
            for j, cursor_val in enumerate(gene.cursor):
                html_file.write("<td>&nbsp;&nbsp;</td>")
                if cursor_val == 0:
                    html_file.write("<td>&nbsp;</td>")
                elif cursor_val == 1:
                    html_file.write("<td>|&nbsp;|</td>")
                elif cursor_val > 1:
                    match_idx = cursor_val - 2
                    if match_idx < len(self.context.match_list):
                        gene2 = self.context.match_list[match_idx].gene2
                        html_file.write(f"<td bgcolor='#ffff99'>{gene2}</td>")
                    else:
                        html_file.write("<td>&nbsp;</td>")
                else:
                    match_idx = -cursor_val - 2
                    if match_idx < len(self.context.match_list):
                        gene1 = self.context.match_list[match_idx].gene1
                        html_file.write(f"<td bgcolor='#ffff99'>{gene1}</td>")
                    else:
                        html_file.write("<td>&nbsp;</td>")

            # 填充剩余列
            for j in range(len(gene.cursor), self.max_level):
                html_file.write("<td>&nbsp;</td>")

            html_file.write("</tr>\n")

        if html_file:
            html_file.write("</table></html>")
            html_file.close()

    def msa_main(self, prefix_fn: str):
        """MSA主函数"""
        self.max_level = 1

        # 获取端点并遍历
        endpoints = self.get_endpoints()
        self.traverse(endpoints)

        # 标记串联重复
        self.mark_tandem(prefix_fn)

        # 输出HTML
        print("Writing multiple syntenic blocks to HTML files")
        self.print_html(prefix_fn)