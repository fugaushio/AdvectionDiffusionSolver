import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse

# コマンドライン引数の設定
parser = argparse.ArgumentParser(description="Visualize time evolution from a CSV file.")
parser.add_argument("csv_file", type=str, help="Path to the CSV file")
args = parser.parse_args()

# CSVファイルの読み込み（ヘッダーなし）
csv_file = args.csv_file
data = pd.read_csv(csv_file, header=None)

# データの形状
num_rows, num_cols = data.shape

# x[m] と time[s] を等間隔で生成
x = np.linspace(0, 1, num_cols)  # 列方向: 0 から 1[m]
time = np.linspace(0, 1, num_rows)  # 行方向: 0 から 1[s]

# アニメーション作成
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)

ax.set_xlim(x.min(), x.max())
ax.set_ylim(data.min().min(), data.max().max())
ax.set_xlabel("x [m]")
ax.set_ylabel("Value")
ax.set_title("Time Evolution")

def init():
    line.set_data([], [])
    return line,

def update(frame):
    line.set_data(x, data.iloc[frame].values)
    ax.set_title(f"Time: {time[frame]:.2f} s")  # 現在時刻を表示
    return line,

# 各フレームの間隔を実際の時間に基づいて計算
frame_interval = (1 / num_rows) * 1000  # ミリ秒単位に変換

ani = FuncAnimation(fig, update, frames=len(time), init_func=init, blit=True, interval=frame_interval)

# アニメーションを保存または表示
# ani.save("animation.mp4", fps=10)  # 保存する場合
plt.show()  # 表示する場合
