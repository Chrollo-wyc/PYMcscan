print("Starting main")
import sys
import os

# 仅在 sys.stderr 存在时启用 faulthandler（避免打包后错误）
if sys.stderr is not None:
    import faulthandler
    faulthandler.enable()

# 初始化 COM（多线程公寓）
try:
    import pythoncom
    pythoncom.CoInitializeEx(pythoncom.COINIT_MULTITHREADED)
    print("COM initialized (MTA) with pythoncom")
except ImportError:
    import ctypes
    ctypes.windll.ole32.CoInitializeEx(None, 0)  # 0 = COINIT_MULTITHREADED
    print("COM initialized (MTA) with ctypes")

# 设置 Qt 插件路径（根据您之前找到的路径）
def get_qt_plugin_path():
    # 根据您实际环境调整（确保路径正确）
    candidate1 = r"E:\anaconda3\envs\pycharm\Library\plugins"
    candidate2 = r"E:\anaconda3\envs\pycharm\Library\lib\qt6\plugins"
    for cand in [candidate1, candidate2]:
        if os.path.exists(cand) and os.path.isdir(cand):
            platforms_dir = os.path.join(cand, 'platforms')
            if os.path.exists(platforms_dir):
                parent = os.path.dirname(cand)
                print(f"找到插件目录: {cand}")
                return parent
    # 备选：PySide6 标准位置
    import site
    for path in site.getsitepackages():
        candidate = os.path.join(path, 'PySide6', 'plugins')
        if os.path.exists(candidate) and os.path.isdir(candidate):
            return path
    return None

plugin_parent = get_qt_plugin_path()
if plugin_parent:
    os.environ['QT_PLUGIN_PATH'] = plugin_parent
    print(f"QT_PLUGIN_PATH set to {plugin_parent}")
else:
    print("错误：未找到 Qt 插件路径，请手动指定")
    sys.exit(1)

# 项目路径设置
from utils import resource_path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from PySide6.QtWidgets import QApplication
from main_window import MainWindow

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()