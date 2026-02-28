#!/usr/bin/env python3
"""
辅助脚本：用于在打包后通过子进程调用 MCScanX 主程序。
接收所有命令行参数，并传递给 pymcscan.cli.mcscan.main()。
"""
import sys
import os

# 确保可以导入本地包
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pymcscan.cli.mcscan import main

if __name__ == "__main__":
    # 注意：main() 函数通常依赖于 sys.argv，因此我们直接传递命令行参数
    # 但也可以选择将参数列表传递给 main 函数（如果它接受参数）
    # 这里假设 main() 使用 sys.argv，因此我们只需调用它
    sys.exit(main())