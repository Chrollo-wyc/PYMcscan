# [file name]: mcscanx/utils/file_utils.py
"""
文件操作和通用工具函数
包含错误处理、进度报告、文件操作等通用功能。
"""

import sys
import os
import time
from datetime import datetime
from typing import Optional, TextIO
import logging

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def progress(format_str: str, *args) -> None:
    """打印进度消息"""
    if args:
        message = format_str % args
    else:
        message = format_str
    print(message, file=sys.stdout)


def err(format_str: str, *args) -> None:
    """打印错误消息但不退出"""
    if args:
        message = format_str % args
    else:
        message = format_str
    logger.error(message)


def warn(format_str: str, *args) -> None:
    """打印警告消息"""
    if args:
        message = format_str % args
    else:
        message = format_str
    logger.warning(message)


def errAbort(format_str: str, *args) -> None:
    """打印错误消息并退出程序"""
    if args:
        message = format_str % args
    else:
        message = format_str
    logger.error(message)
    print(f"[Error] {message}", file=sys.stderr)
    sys.exit(1)


class Timer:
    """计时器类"""

    def __init__(self):
        self.start_time = time.time()
        self.last_time = self.start_time

    def elapsed(self) -> float:
        """返回从开始到现在的时间（秒）"""
        return time.time() - self.start_time

    def lap(self) -> float:
        """返回从上次调用lap到现在的时间（秒）"""
        current = time.time()
        elapsed = current - self.last_time
        self.last_time = current
        return elapsed


# 全局计时器实例
_timer = Timer()


def uglyTime(label: Optional[str] = None, *args) -> None:
    """
    打印标签和自上次调用以来的时间
    使用None作为标签进行初始化
    """
    if label is None:
        # 重新初始化计时器
        global _timer
        _timer = Timer()
        return

    if args:
        formatted_label = label % args
    else:
        formatted_label = label

    elapsed = _timer.lap()
    print(f"{formatted_label} [{elapsed:.3f} seconds elapsed]")


def clock1000() -> int:
    """返回毫秒级时钟"""
    return int(time.time() * 1000)


def mustOpen(fileName: str, mode: str = 'r') -> TextIO:
    """
    打开文件，如果失败则退出程序
    支持特殊的文件名 'stdin', 'stdout', 'stderr'
    """
    if fileName == "stdin":
        return sys.stdin
    if fileName == "stdout":
        return sys.stdout
    if fileName == "stderr":
        return sys.stderr

    try:
        if 'w' in mode or 'a' in mode:
            # 对于写模式，确保目录存在
            dirname = os.path.dirname(fileName)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname, exist_ok=True)

        return open(fileName, mode)
    except IOError as e:
        mode_name = ""
        if mode:
            if mode[0] == 'r':
                mode_name = " to read"
            elif mode[0] == 'w':
                mode_name = " to write"
            elif mode[0] == 'a':
                mode_name = " to append"

        errAbort("Can't open %s%s: %s", fileName, mode_name, str(e))


def sameString(a: str, b: str) -> bool:
    """比较两个字符串是否相同"""
    return a == b


def file_exists(filepath: str) -> bool:
    """检查文件是否存在且可读"""
    return os.path.exists(filepath) and os.access(filepath, os.R_OK)


def create_directory(dirpath: str) -> bool:
    """创建目录，如果已存在则返回True"""
    try:
        os.makedirs(dirpath, exist_ok=True)
        return True
    except OSError as e:
        err(f"Failed to create directory {dirpath}: {str(e)}")
        return False


def get_file_size(filepath: str) -> int:
    """获取文件大小（字节）"""
    try:
        return os.path.getsize(filepath)
    except OSError:
        return 0


def read_lines(filepath: str, skip_comments: bool = True, comment_char: str = '#') -> list:
    """读取文件行，可选跳过注释行"""
    lines = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if skip_comments and line.startswith(comment_char):
                    continue
                if line:  # 跳过空行
                    lines.append(line)
    except IOError as e:
        err(f"Error reading file {filepath}: {str(e)}")
    return lines


def write_lines(filepath: str, lines: list, header: str = None) -> bool:
    """写入多行到文件"""
    try:
        with open(filepath, 'w') as f:
            if header:
                f.write(header + '\n')
            for line in lines:
                f.write(str(line) + '\n')
        return True
    except IOError as e:
        err(f"Error writing to file {filepath}: {str(e)}")
        return False


def format_number(n: int) -> str:
    """格式化大数字，添加千位分隔符"""
    return f"{n:,}"


def format_percentage(part: int, total: int) -> str:
    """格式化百分比"""
    if total == 0:
        return "0.00%"
    return f"{100.0 * part / total:.2f}%"


def get_timestamp() -> str:
    """获取当前时间戳字符串"""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class ProgressBar:
    """进度条类"""

    def __init__(self, total: int, width: int = 50, desc: str = "Progress"):
        self.total = total
        self.width = width
        self.desc = desc
        self.current = 0

    def update(self, increment: int = 1) -> None:
        """更新进度"""
        self.current += increment
        if self.current > self.total:
            self.current = self.total

        percent = self.current / self.total
        filled = int(self.width * percent)
        bar = '█' * filled + '░' * (self.width - filled)

        sys.stdout.write(f"\r{self.desc}: |{bar}| {percent:.1%} ({self.current}/{self.total})")
        sys.stdout.flush()

    def finish(self) -> None:
        """完成进度条"""
        self.update(self.total - self.current)
        print()  # 换行