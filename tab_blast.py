import os
import subprocess
from input import ParamInput, FileGroupBox
from PySide6.QtWidgets import (
    QMessageBox, QProgressDialog, QPushButton,
    QWidget, QVBoxLayout, QGroupBox, QGridLayout, QLabel,
    QScrollArea
)
from PySide6.QtCore import QThread, Qt, Signal
from utils import resource_path


class BlastParamGroupBox(QGroupBox):
    """BLAST 参数设置盒子，带标签的水平网格布局"""
    def __init__(self, title="BLAST 参数设置", parent=None):
        super().__init__(title, parent)
        layout = QGridLayout()

        # 第一行：BLAST类型、Query类型、Database类型
        layout.addWidget(QLabel("BLAST类型:"), 0, 0)
        self.param1 = ParamInput("blastp")
        self.param1.setText("blastp")
        layout.addWidget(self.param1, 0, 1)

        layout.addWidget(QLabel("Query类型:"), 0, 2)
        self.param2 = ParamInput("prot")
        self.param2.setText("prot")
        layout.addWidget(self.param2, 0, 3)

        layout.addWidget(QLabel("Database类型:"), 0, 4)
        self.param3 = ParamInput("prot")
        self.param3.setText("prot")
        layout.addWidget(self.param3, 0, 5)

        # 第二行：E值、输出格式、线程数
        layout.addWidget(QLabel("E-value阈值:"), 1, 0)
        self.param4 = ParamInput("1e-5")
        self.param4.setText("1e-5")
        layout.addWidget(self.param4, 1, 1)

        layout.addWidget(QLabel("输出格式:"), 1, 2)
        self.param5 = ParamInput("8")
        self.param5.setText("8")
        layout.addWidget(self.param5, 1, 3)

        layout.addWidget(QLabel("线程数:"), 1, 4)
        self.param6 = ParamInput("4")
        self.param6.setText("4")
        layout.addWidget(self.param6, 1, 5)

        # 第三行：覆盖度
        layout.addWidget(QLabel("覆盖度:"), 2, 0)
        self.param7 = ParamInput("80%")
        self.param7.setText("80%")
        layout.addWidget(self.param7, 2, 1)

        # 让输入框所在列自动拉伸
        layout.setColumnStretch(1, 1)
        layout.setColumnStretch(3, 1)
        layout.setColumnStretch(5, 1)

        self.setLayout(layout)


class BuildDBThread(QThread):
    """建库工作线程"""
    finished = Signal(bool, str)  # 成功标志，消息
    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd

    def run(self):
        try:
            subprocess.run(self.cmd, check=True, capture_output=True, text=True, encoding='utf-8')
            self.finished.emit(True, "建库成功")
        except subprocess.CalledProcessError as e:
            self.finished.emit(False, f"建库失败:\n{e.stderr}")
        except Exception as e:
            self.finished.emit(False, f"建库异常:\n{str(e)}")


class BlastRunThread(QThread):
    """BLAST比对工作线程"""
    finished = Signal(bool, str)
    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd

    def run(self):
        try:
            subprocess.run(self.cmd, check=True, capture_output=True, text=True, encoding='utf-8')
            self.finished.emit(True, "比对完成")
        except subprocess.CalledProcessError as e:
            self.finished.emit(False, f"比对失败:\n{e.stderr}")
        except Exception as e:
            self.finished.emit(False, f"比对异常:\n{str(e)}")


class BlastTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.build_thread = None
        self.blast_thread = None
        self.progress = None
        self.initUI()

    def initUI(self):
        # 创建内容容器布局
        layout = QVBoxLayout()

        # 文件选择盒子
        self.file1 = FileGroupBox("Query sequence")
        self.file2 = FileGroupBox("Database sequence")
        self.db_path = FileGroupBox("Database path")
        self.result_path = FileGroupBox("Result path")

        # 参数盒子
        self.blast_param = BlastParamGroupBox()

        # 开始按钮
        self.start_button = QPushButton("START")
        self.start_button.clicked.connect(self.run_blast)

        # 将所有部件添加到布局
        for widget in [self.file1, self.file2, self.db_path,
                       self.result_path, self.blast_param, self.start_button]:
            layout.addWidget(widget)

        layout.addStretch()

        # 创建滚动区域，将内容放入其中
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        content_widget = QWidget()
        content_widget.setLayout(layout)
        scroll_area.setWidget(content_widget)

        # 将滚动区域作为标签页的主布局
        main_layout = QVBoxLayout()
        main_layout.addWidget(scroll_area)
        self.setLayout(main_layout)

    def get_blast_exe(self, blast_type):
        exe_name = f"{blast_type}.exe" if os.name == 'nt' else blast_type
        relative_path = os.path.join("blast", "bin", exe_name)
        return resource_path(relative_path)

    # 注意：原来的 create_blast_db 方法已不再使用，可删除或注释掉
    # def create_blast_db(self, db, db_path, file2_type):
    #     ...

    def run_blast(self):
        # ---------- 获取并验证输入 ----------
        db_path = self.db_path.path_input.text().strip()
        result_path = self.result_path.path_input.text().strip()
        query = self.file1.path_input.text().strip()
        db_seq = self.file2.path_input.text().strip()   # 数据库序列文件
        blast_type = self.blast_param.param1.text().strip()

        if not db_path or not result_path or not query or not db_seq or not blast_type:
            QMessageBox.warning(self, "警告", "请填写所有路径和参数")
            return

        if not os.path.isfile(query):
            QMessageBox.critical(self, "错误", f"查询文件不存在：{query}")
            return
        if not os.path.isfile(db_seq):
            QMessageBox.critical(self, "错误", f"数据库序列文件不存在：{db_seq}")
            return

        # 确保输出目录存在
        try:
            os.makedirs(db_path, exist_ok=True)
            os.makedirs(result_path, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(self, "错误", f"无法创建输出目录：{str(e)}")
            return

        # 定义数据库索引文件列表
        db_index_files = ['db.pin', 'db.phr', 'db.psq', 'db.pot']
        db_exists = all(os.path.exists(os.path.join(db_path, f)) for f in db_index_files)

        # 如果数据库不存在，询问是否建库
        if not db_exists:
            reply = QMessageBox.question(
                self, "数据库不存在",
                f"数据库文件在 {db_path} 中不存在，是否根据序列文件 {db_seq} 创建？",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.No:
                return
            # 启动建库线程
            self.start_build_db(db_seq, db_path)
        else:
            # 数据库已存在，直接启动比对
            self.start_blast_run(query, db_path, result_path, blast_type)

    def start_build_db(self, db_seq, db_path):
        """启动建库线程"""
        # 获取 makeblastdb 命令
        makeblastdb_exe = self.get_blast_exe("makeblastdb")
        cmd = [
            makeblastdb_exe,
            '-in', db_seq,
            '-out', os.path.join(db_path, 'db'),
            '-dbtype', self.blast_param.param3.text().strip(),
            '-parse_seqids'
        ]

        # 显示无限进度对话框
        self.progress = QProgressDialog("正在创建BLAST数据库...", "取消", 0, 0, self)
        self.progress.setWindowTitle("进度")
        self.progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress.setMinimumDuration(0)
        self.progress.canceled.connect(self.cancel_build_db)

        self.build_thread = BuildDBThread(cmd)
        self.build_thread.finished.connect(self.on_build_finished)
        self.build_thread.start()
        self.progress.show()

    def cancel_build_db(self):
        if hasattr(self, 'build_thread') and self.build_thread.isRunning():
            self.build_thread.terminate()
            self.build_thread.wait()
        self.progress.close()
        QMessageBox.information(self, "取消", "建库已取消")

    def on_build_finished(self, success, message):
        self.progress.close()
        if success:
            # 建库成功，启动比对
            self.start_blast_run(
                self.file1.path_input.text().strip(),
                self.db_path.path_input.text().strip(),
                self.result_path.path_input.text().strip(),
                self.blast_param.param1.text().strip()
            )
        else:
            QMessageBox.critical(self, "建库失败", message)
        self.build_thread = None

    def start_blast_run(self, query, db_path, result_path, blast_type):
        """启动比对线程"""
        output_file = os.path.join(result_path, "results.blast")
        db_file = os.path.join(db_path, 'db')

        blast_exe = self.get_blast_exe(blast_type)
        cmd = [
            blast_exe,
            '-query', query,
            '-db', db_file,
            '-out', output_file,
            '-evalue', self.blast_param.param4.text().strip(),
            '-outfmt', self.blast_param.param5.text().strip(),
            '-num_threads', self.blast_param.param6.text().strip()
        ]

        self.progress = QProgressDialog("正在执行BLAST比对...", "取消", 0, 0, self)
        self.progress.setWindowTitle("进度")
        self.progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress.setMinimumDuration(0)
        self.progress.canceled.connect(self.cancel_blast_run)

        self.blast_thread = BlastRunThread(cmd)
        self.blast_thread.finished.connect(self.on_blast_finished)
        self.blast_thread.start()
        self.progress.show()

    def cancel_blast_run(self):
        if hasattr(self, 'blast_thread') and self.blast_thread.isRunning():
            self.blast_thread.terminate()
            self.blast_thread.wait()
        self.progress.close()
        QMessageBox.information(self, "取消", "比对已取消")

    def on_blast_finished(self, success, message):
        self.progress.close()
        if success:
            QMessageBox.information(self, "成功", f"BLAST分析完成！\n结果保存至: {self.result_path.path_input.text()}")
        else:
            QMessageBox.critical(self, "比对失败", message)
        self.blast_thread = None