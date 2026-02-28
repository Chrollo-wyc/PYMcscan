import os
import sys
import traceback
from input import ParamInput, FileGroupBox
from PySide6.QtWidgets import (
    QMessageBox, QProgressDialog, QPushButton,
    QWidget, QVBoxLayout, QGroupBox, QGridLayout, QLabel,
    QScrollArea, QCheckBox
)
from PySide6.QtCore import QThread, Qt, Signal
from mycircos.draw import draw_circos


# PyCircos 参数设置盒子
class PyCircosParamGroupBox(QGroupBox):
    def __init__(self, title="PyCircos 参数设置", parent=None):
        super().__init__(title, parent)
        layout = QGridLayout()

        # 第一行：图标题、宽度、高度
        layout.addWidget(QLabel("图标题:"), 0, 0)
        self.param_title = ParamInput("Genome Circos Plot", "Genome Circos Plot")
        layout.addWidget(self.param_title, 0, 1)

        layout.addWidget(QLabel("宽度(英寸):"), 0, 2)
        self.param_figwidth = ParamInput("10", "10")
        layout.addWidget(self.param_figwidth, 0, 3)

        layout.addWidget(QLabel("高度(英寸):"), 0, 4)
        self.param_figheight = ParamInput("8", "8")
        layout.addWidget(self.param_figheight, 0, 5)

        # 第二行：DPI、连接线颜色、透明度
        layout.addWidget(QLabel("DPI:"), 1, 0)
        self.param_dpi = ParamInput("300", "300")
        layout.addWidget(self.param_dpi, 1, 1)

        layout.addWidget(QLabel("连接线颜色:"), 1, 2)
        self.param_link_color = ParamInput("red", "red")
        layout.addWidget(self.param_link_color, 1, 3)

        layout.addWidget(QLabel("透明度:"), 1, 4)
        self.param_link_alpha = ParamInput("0.3", "0.3")
        layout.addWidget(self.param_link_alpha, 1, 5)

        # 第三行：是否显示标签（可选）
        self.param_show_labels = QCheckBox("显示染色体标签")
        self.param_show_labels.setChecked(True)
        layout.addWidget(self.param_show_labels, 2, 0, 1, 2)

        # 让输入框所在列自动拉伸
        for col in [1, 3, 5]:
            layout.setColumnStretch(col, 1)

        self.setLayout(layout)


# 工作线程：直接导入 draw 模块并调用绘图函数
class PyCircosWorker(QThread):
    finished = Signal(int, str)  # 返回码和输出信息

    def __init__(self, params):
        super().__init__()
        self.params = params

    def run(self):
        try:
            # 导入绘图模块（确保 mycircos.draw 在 Python 路径中）


            # 提取参数（并转换类型）
            gff_file = self.params['gff_file']
            coll_file = self.params['coll_file']
            output_file = self.params['output_file']
            title = self.params['title']
            figwidth = float(self.params['figwidth'])
            figheight = float(self.params['figheight'])
            dpi = int(self.params['dpi'])
            link_color = self.params['link_color']
            link_alpha = float(self.params['link_alpha'])
            show_labels = self.params.get('show_labels', True)

            # 调用绘图函数（请确保函数签名与此一致）
            draw_circos(
                gff_file=gff_file,
                collinearity_file=coll_file,
                output_file=output_file,
                title=title,
                figsize=(figwidth, figheight),
                dpi=dpi,
                link_color=link_color,
                link_alpha=link_alpha,
                show_labels=show_labels
            )

            self.finished.emit(0, "绘图成功")

        except ImportError as e:
            self.finished.emit(-1, f"无法导入 mycircos.draw 模块: {str(e)}\n请确保 mycircos 包已正确安装，且 draw.py 文件存在。")
        except Exception as e:
            error_msg = f"绘图失败: {str(e)}\n{traceback.format_exc()}"
            self.finished.emit(-1, error_msg)


# PyCircos 标签页主类
class PyCircosTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        # 文件选择
        self.gff_box = FileGroupBox("GFF file (.gff)")
        layout.addWidget(self.gff_box)

        self.collinearity_box = FileGroupBox("Collinearity file (.collinearity)")
        layout.addWidget(self.collinearity_box)

        self.output_dir_box = FileGroupBox("Output directory")
        layout.addWidget(self.output_dir_box)

        self.prefix_input = ParamInput("Output file prefix (e.g., circos_plot)", "circos_plot")
        layout.addWidget(self.prefix_input)

        # 参数设置
        self.param_box = PyCircosParamGroupBox()
        layout.addWidget(self.param_box)

        # 运行按钮
        self.run_button = QPushButton("RUN PyCircos")
        self.run_button.clicked.connect(self.run_pycircos)
        layout.addWidget(self.run_button)

        layout.addStretch()

        # 滚动区域
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        content_widget = QWidget()
        content_widget.setLayout(layout)
        scroll_area.setWidget(content_widget)

        main_layout = QVBoxLayout()
        main_layout.addWidget(scroll_area)
        self.setLayout(main_layout)

        self.worker = None
        self.progress = None

    def run_pycircos(self):
        # 获取用户输入
        gff_file = self.gff_box.path_input.text().strip()
        coll_file = self.collinearity_box.path_input.text().strip()
        output_dir = self.output_dir_box.path_input.text().strip()
        prefix = self.prefix_input.text().strip()

        if not gff_file or not coll_file or not output_dir or not prefix:
            QMessageBox.warning(self, "警告", "请填写所有文件路径、输出目录和前缀")
            return

        if not os.path.isfile(gff_file):
            QMessageBox.critical(self, "错误", f"GFF 文件不存在：{gff_file}")
            return
        if not os.path.isfile(coll_file):
            QMessageBox.critical(self, "错误", f"共线性文件不存在：{coll_file}")
            return

        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(self, "错误", f"无法创建输出目录：{str(e)}")
            return

        # 输出文件路径
        output_file = os.path.join(output_dir, f"{prefix}.png")

        # 收集参数
        title = self.param_box.param_title.text().strip()
        figwidth = self.param_box.param_figwidth.text().strip()
        figheight = self.param_box.param_figheight.text().strip()
        dpi = self.param_box.param_dpi.text().strip()
        link_color = self.param_box.param_link_color.text().strip()
        link_alpha = self.param_box.param_link_alpha.text().strip()
        show_labels = self.param_box.param_show_labels.isChecked()

        params = {
            'gff_file': gff_file,
            'coll_file': coll_file,
            'output_file': output_file,
            'title': title,
            'figwidth': figwidth,
            'figheight': figheight,
            'dpi': dpi,
            'link_color': link_color,
            'link_alpha': link_alpha,
            'show_labels': show_labels
        }

        self.progress = QProgressDialog("正在绘制 Circos 图...", "取消", 0, 0, self)
        self.progress.setWindowTitle("进度")
        self.progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress.setMinimumDuration(0)
        self.progress.canceled.connect(self.cancel_worker)

        self.worker = PyCircosWorker(params)
        self.worker.finished.connect(self.on_worker_finished)
        self.worker.start()
        self.progress.show()

    def cancel_worker(self):
        if self.worker and self.worker.isRunning():
            self.worker.terminate()
            self.worker.wait()
        self.progress.close()
        QMessageBox.information(self, "取消", "绘图已取消")

    def on_worker_finished(self, returncode, output):
        self.progress.close()
        if returncode == 0:
            QMessageBox.information(
                self, "成功",
                f"Circos 图绘制完成！\n结果文件已生成于：{self.output_dir_box.path_input.text()}"
            )
        else:
            short_output = output[:500] + ("..." if len(output) > 500 else "")
            QMessageBox.critical(
                self, "错误",
                f"绘图失败，返回码 {returncode}\n错误输出：\n{short_output}"
            )
        self.worker = None