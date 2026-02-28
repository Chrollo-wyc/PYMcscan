import os
import sys
import shutil
import subprocess
from input import ParamInput, FileGroupBox
from PySide6.QtWidgets import (
    QMessageBox, QProgressDialog, QPushButton,
    QWidget, QVBoxLayout, QGroupBox, QGridLayout, QLabel,
    QCheckBox, QScrollArea
)
from PySide6.QtCore import QThread, Qt, Signal
from utils import resource_path


# MCScanX 参数设置盒子
class MCScanXParamGroupBox(QGroupBox):
    def __init__(self, title="MCScanX 参数设置", parent=None):
        super().__init__(title, parent)
        layout = QGridLayout()

        # 第一行：-k, -g, -s
        layout.addWidget(QLabel("Match score (-k):"), 0, 0)
        self.param_k = ParamInput("50")
        self.param_k.setText("50")
        layout.addWidget(self.param_k, 0, 1)

        layout.addWidget(QLabel("Gap penalty (-g):"), 0, 2)
        self.param_g = ParamInput("-1")
        self.param_g.setText("-1")
        layout.addWidget(self.param_g, 0, 3)

        layout.addWidget(QLabel("Match size (-s):"), 0, 4)
        self.param_s = ParamInput("5")
        self.param_s.setText("5")
        layout.addWidget(self.param_s, 0, 5)

        # 第二行：-e, -m, -w
        layout.addWidget(QLabel("E-value (-e):"), 1, 0)
        self.param_e = ParamInput("1e-5")
        self.param_e.setText("1e-5")
        layout.addWidget(self.param_e, 1, 1)

        layout.addWidget(QLabel("Max gaps (-m):"), 1, 2)
        self.param_m = ParamInput("25")
        self.param_m.setText("25")
        layout.addWidget(self.param_m, 1, 3)

        layout.addWidget(QLabel("Overlap window (-w):"), 1, 4)
        self.param_w = ParamInput("5")
        self.param_w.setText("5")
        layout.addWidget(self.param_w, 1, 5)

        # 第三行：-b, -c, 以及复选框
        layout.addWidget(QLabel("Block pattern (-b):"), 2, 0)
        self.param_b = ParamInput("0")
        self.param_b.setText("0")
        layout.addWidget(self.param_b, 2, 1)

        layout.addWidget(QLabel("Homology score (-c):"), 2, 2)
        self.param_c = ParamInput("0")
        self.param_c.setText("0")
        layout.addWidget(self.param_c, 2, 3)

        self.check_a = QCheckBox("Pairwise only (-a)")
        layout.addWidget(self.check_a, 2, 4, 1, 2)

        for col in [1, 3, 5]:
            layout.setColumnStretch(col, 1)
        self.setLayout(layout)


# MCScanX 工作线程
class MCScanXWorker(QThread):
    finished = Signal(int, str)

    def __init__(self, args):
        """
        args: 命令行参数列表（不包括脚本名），例如 ['prefix', '-k', '50', ...]
        """
        super().__init__()
        self.args = args

    def run(self):
        import traceback
        log_file = os.path.join(os.path.expanduser("~"), "mcscan_worker.log")
        with open(log_file, "a") as log:
            log.write(f"Worker started with args: {self.args}\n")
            original_argv = sys.argv
            try:
                sys.argv = [sys.argv[0]] + self.args
                log.write(f"sys.argv set to: {sys.argv}\n")
                from pymcscan.cli.mcscan import main
                log.write("Imported main successfully\n")
                try:
                    main()
                    returncode = 0
                    output = ""
                    log.write("main() completed normally\n")
                except SystemExit as e:
                    returncode = e.code if isinstance(e.code, int) else 1
                    output = str(e) if str(e) else f"MCScanX exited with code {returncode}"
                    log.write(f"SystemExit caught: {returncode}\n")
                except Exception as e:
                    returncode = -1
                    output = f"Unexpected error: {str(e)}\n{traceback.format_exc()}"
                    log.write(f"Exception in main: {str(e)}\n{traceback.format_exc()}\n")
            except ImportError as e:
                returncode = -1
                output = f"Failed to import pymcscan: {str(e)}\n{traceback.format_exc()}"
                log.write(f"ImportError: {output}\n")
            except Exception as e:
                returncode = -1
                output = f"Critical error in worker: {str(e)}\n{traceback.format_exc()}"
                log.write(f"Critical error: {output}\n")
            finally:
                sys.argv = original_argv
                log.write("Restored sys.argv\n")
            log.write(f"Emitting finished signal with returncode {returncode}\n")
            self.finished.emit(returncode, output)


# MCScanX 标签页主类
class MCScanXTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()

    def initUI(self):
        # 内容布局
        layout = QVBoxLayout()

        # 文件选择
        self.gff_box = FileGroupBox("GFF file (.gff)")
        layout.addWidget(self.gff_box)

        self.blast_box = FileGroupBox("BLAST file (.blast or .homology)")
        layout.addWidget(self.blast_box)

        self.output_dir_box = FileGroupBox("Output directory")
        layout.addWidget(self.output_dir_box)

        self.prefix_input = ParamInput("Output file prefix (e.g., mcscan_result)", "mcscan_result")
        layout.addWidget(self.prefix_input)

        # 参数设置
        self.param_box = MCScanXParamGroupBox()
        layout.addWidget(self.param_box)

        # 运行按钮
        self.run_button = QPushButton("RUN MCScanX")
        self.run_button.clicked.connect(self.run_mcscan)
        layout.addWidget(self.run_button)

        layout.addStretch()

        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        content_widget = QWidget()
        content_widget.setLayout(layout)
        scroll_area.setWidget(content_widget)

        # 主布局
        main_layout = QVBoxLayout()
        main_layout.addWidget(scroll_area)
        self.setLayout(main_layout)

        self.mcscan_worker = None
        self.progress = None

    def run_mcscan(self):
        gff_src = self.gff_box.path_input.text().strip()
        blast_src = self.blast_box.path_input.text().strip()
        output_dir = self.output_dir_box.path_input.text().strip()
        prefix = self.prefix_input.text().strip()

        if not gff_src or not blast_src or not output_dir or not prefix:
            QMessageBox.warning(self, "警告", "请填写所有文件路径、输出目录和前缀")
            return

        if not os.path.isfile(gff_src):
            QMessageBox.critical(self, "错误", f"GFF 文件不存在：{gff_src}")
            return
        if not os.path.isfile(blast_src):
            QMessageBox.critical(self, "错误", f"BLAST 文件不存在：{blast_src}")
            return

        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(self, "错误", f"无法创建输出目录：{str(e)}")
            return

        target_gff = os.path.join(output_dir, prefix + ".gff")
        target_blast = os.path.join(output_dir, prefix + ".blast")

        def safe_copy(src, dst, desc):
            src_abs = os.path.abspath(src)
            dst_abs = os.path.abspath(dst)
            if src_abs == dst_abs:
                return True
            if os.path.exists(dst_abs):
                reply = QMessageBox.question(
                    self, "文件已存在",
                    f"{desc} 目标文件已存在：\n{dst_abs}\n是否覆盖？",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply == QMessageBox.StandardButton.No:
                    return False
            try:
                shutil.copy2(src, dst)
                return True
            except Exception as e:
                QMessageBox.critical(self, "错误", f"复制 {desc} 失败：{str(e)}")
                return False

        if not safe_copy(gff_src, target_gff, "GFF 文件"):
            return
        if not safe_copy(blast_src, target_blast, "BLAST 文件"):
            return

        full_prefix = os.path.join(output_dir, prefix)

        # 收集 MCScanX 参数（不包括脚本名）
        args = [full_prefix]

        if self.param_box.param_k.text():
            args += ["-k", self.param_box.param_k.text()]
        if self.param_box.param_g.text():
            args += ["-g", self.param_box.param_g.text()]
        if self.param_box.param_s.text():
            args += ["-s", self.param_box.param_s.text()]
        if self.param_box.param_e.text():
            args += ["-e", self.param_box.param_e.text()]
        if self.param_box.param_m.text():
            args += ["-m", self.param_box.param_m.text()]
        if self.param_box.param_w.text():
            args += ["-w", self.param_box.param_w.text()]
        if self.param_box.param_b.text():
            args += ["-b", self.param_box.param_b.text()]
        if self.param_box.param_c.text():
            args += ["-c", self.param_box.param_c.text()]
        if self.param_box.check_a.isChecked():
            args.append("-a")

        self.progress = QProgressDialog("正在运行 MCScanX...", "取消", 0, 0, self)
        self.progress.setWindowTitle("进度")
        self.progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress.setMinimumDuration(0)
        self.progress.canceled.connect(self.cancel_mcscan)

        # 创建并启动工作线程（传递 args）
        self.mcscan_worker = MCScanXWorker(args)
        self.mcscan_worker.finished.connect(self.on_mcscan_finished)
        self.mcscan_worker.start()
        self.progress.show()

    def cancel_mcscan(self):
        if self.mcscan_worker and self.mcscan_worker.isRunning():
            self.mcscan_worker.terminate()
            self.mcscan_worker.wait()
        self.progress.close()
        QMessageBox.information(self, "取消", "MCScanX 运行已取消")

    def on_mcscan_finished(self, returncode, output):
        log_file = os.path.join(os.path.expanduser("~"), "mcscan_worker.log")
        with open(log_file, "a") as log:
            log.write(f"on_mcscan_finished called with returncode {returncode}, output length: {len(output)}\n")
        self.progress.close()
        if returncode == 0:
            QMessageBox.information(
                self, "成功",
                f"MCScanX 分析完成！\n结果文件已生成于：{self.output_dir_box.path_input.text()}"
            )
        else:
            short_output = output[:500] + ("..." if len(output) > 500 else "")
            QMessageBox.critical(
                self, "错误",
                f"MCScanX 运行失败，返回码 {returncode}\n错误输出：\n{short_output}"
            )
        self.mcscan_worker = None
