from PySide6.QtWidgets import QMainWindow, QTabWidget
from tab_blast import BlastTab
from tab_pymcscan import MCScanXTab
from tab_pycircos import PyCircosTab

#主窗口
class MainWindow(QMainWindow):
    #初始化
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('PyMcscan')
        self.setGeometry(100, 100, 800, 600)
        # 创建“标签页容器”
        self.tab_widget = QTabWidget(self)
        self.setCentralWidget(self.tab_widget)
    # 添加标签页
        # 创建“blast”标签页
        tab_blast = BlastTab()
        self.tab_widget.addTab(tab_blast, "Blast")
        # 创建“mcscan”标签页
        tab_pymcscan =MCScanXTab()
        self.tab_widget.addTab(tab_pymcscan, "MCScanX")
        # 创建“CIRCOS”标签页
        tab_pycircos = PyCircosTab()
        self.tab_widget.addTab(tab_pycircos, "PyCircos")
