import os
import subprocess


def create_blast_db(fasta_path, db_path, db_type='prot'):
    """
    使用makeblastdb创建BLAST数据库
    :param fasta_path: 输入FASTA文件路径
    :param db_path: 输出数据库路径
    :param db_type: 数据库类型 ('prot' 或 'nucl')
    """
    # 验证输入文件
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"输入文件不存在: {fasta_path}")

        # 创建输出目录
    os.makedirs(db_path, exist_ok=True)

    # 构建命令
    cmd = [
        'D:/zip/D/blast/bin/makeblastdb.exe',
        '-in', fasta_path,
        '-out', os.path.join(db_path, 'db'),  # 数据库文件名
        '-dbtype', db_type,
        '-parse_seqids'
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"BLAST数据库创建成功，路径: {os.path.join(db_path, 'db')}")
    except subprocess.CalledProcessError as e:
        print(f"命令输出:\n{e.stdout}")
        print(f"错误输出:\n{e.stderr}")
        raise RuntimeError("BLAST数据库创建失败")


def run_blastp(fasta_path, db_path, output_path):
    """
    执行BLASTp比对
    :param fasta_path: 查询序列路径
    :param db_path: BLAST数据库路径
    :param output_path: 输出文件路径
    """
    # 验证数据库存在性
    db_files = ['db.pin', 'db.phr', 'db.psq', 'db.pot']
    db_exists = all(os.path.exists(os.path.join(db_path, f)) for f in db_files)
    if not db_exists:
        raise FileNotFoundError(f"BLAST数据库文件在 {db_path} 中缺失")

        # 构建命令
    cmd = [
        'D:/zip/D/blast/bin/blastp.exe',
        '-query', fasta_path,
        '-db', os.path.join(db_path, 'db'),
        '-out', output_path,
        '-evalue', '0.001',
        '-outfmt', '6',
        '-num_threads', '4'  # 启用多线程
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"BLASTp比对完成，结果保存至: {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"命令输出:\n{e.stdout}")
        print(f"错误输出:\n{e.stderr}")
        raise RuntimeError("BLASTp比对失败")



# 主程序
if __name__ == '__main__':
    # 配置参数
    fasta_path = "D:/zip/D/blast/Gma.fa"
    db_fasta = "D:/zip/D/blast/Ath.fa"
    db_path = "D:/zip/D/blast/"
    output_path = "D:/zip/D/blast/results.blast"

    try:
        # 创建数据库
        create_blast_db(db_fasta, db_path, 'prot')

        # 执行比对
        run_blastp(fasta_path, db_path, output_path)
    except Exception as e:
        print(f"处理失败: {str(e)}")
