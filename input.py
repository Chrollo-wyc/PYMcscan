from PySide6.QtWidgets import QFileDialog, QVBoxLayout, QGroupBox,QLabel,QLineEdit,QPushButton


#创建一个’文件选择盒子‘类
class FileGroupBox(QGroupBox):
    def __init__(self, title="文件选择", parent=None):
        super().__init__(title, parent)  # 设置盒子标题
        # 第二步：搭“内部小隔板”（垂直布局）
        self.inner_layout = QVBoxLayout()
        # 第三步：往盒子里装“零件”
        self.title_label = QLabel(f"请选择{title}：")
        self.path_input = QLineEdit()
        self.browse_button = QPushButton("浏览")

        self.inner_layout.addWidget(self.title_label)
        self.inner_layout.addWidget(self.path_input)
        self.inner_layout.addWidget(self.browse_button)

        self.setLayout(self.inner_layout)

        self.browse_button.clicked.connect(self.open_file_dialog)
    # 定义打开文件对话框的函数
    def open_file_dialog(self):
        # 弹出文件选择对话框，返回选择的文件路径
        file_path, _ = QFileDialog.getOpenFileName(self, "选择文件")
        if file_path:  # 如果选择了文件
            self.path_input.setText(file_path)  # 在输入框显示路径


#创建输入框类
class ParamInput(QLineEdit):
    def __init__(self, title="请输入参数：", default_value=""):
        super().__init__()  # 调用爸爸的构造方法
        self.setPlaceholderText(title)  # 设置输入框提示文本（灰色占位符）
        self.setText(default_value)  # 设置默认值
        self.setClearButtonEnabled(True)  # 启用右侧清除按钮（点×清空内容）