
# 导入networkx库
import networkx as nx
import numpy as np
# 定义一个函数，用来生成一个随机网络
def generate_random_network(n, p):
    # n是节点数，p是连边概率
    # 返回一个networkx.Graph对象
    return nx.erdos_renyi_graph(n, p)

# 定义一个函数，用来模拟SIS模型下的疾病传播
def simulate_SIS(G, beta, gamma, T):
    # G是网络，beta是传染率，gamma是恢复率，T是时间步数
    # 返回一个列表，表示每个时间步下的感染节点数
    # 初始化所有节点为易感状态
    nx.set_node_attributes(G, 0, 'state')
    # 随机选择一个节点作为初始感染者
    infected_node = np.random.choice(G.nodes())
    G.nodes[infected_node]['state'] = 1
    # 初始化感染节点数列表
    infected_count = [1]
    # 循环T次
    for t in range(T):
        # 遍历所有节点
        for node in G.nodes():
            # 如果节点是感染状态
            if G.nodes[node]['state'] == 1:
                # 以gamma的概率恢复为易感状态
                if np.random.random() < gamma:
                    G.nodes[node]['state'] = 0
                # 否则，遍历所有邻居节点
                else:
                    for neighbor in G.neighbors(node):
                        # 如果邻居节点是易感状态，以beta的概率感染邻居节点
                        if G.nodes[neighbor]['state'] == 0 and np.random.random() < beta:
                            G.nodes[neighbor]['state'] = 1
        # 计算当前时间步下的感染节点数，并添加到列表中
        infected_count.append(sum(nx.get_node_attributes(G, 'state').values()))
    # 返回感染节点数列表
    return infected_count

# 定义一个函数，用来模拟SIR模型下的疾病传播
def simulate_SIR(G, beta, gamma, T):
    # G是网络，beta是传染率，gamma是恢复率，T是时间步数
    # 返回一个列表，表示每个时间步下的感染节点数
    # 初始化所有节点为易感状态
    nx.set_node_attributes(G, 0, 'state')
    # 随机选择一个节点作为初始感染者
    infected_node = np.random.choice(G.nodes())
    G.nodes[infected_node]['state'] = 1
    # 初始化感染节点数列表
    infected_count = [1]
    # 循环T次
    for t in range(T):
        # 遍历所有节点
        for node in G.nodes():
            # 如果节点是感染状态
            if G.nodes[node]['state'] == 1:
                # 以gamma的概率恢复为免疫状态
                if np.random.random() < gamma:
                    G.nodes[node]['state'] = 2
                # 否则，遍历所有邻居节点
                else:
                    for neighbor in G.neighbors(node):
                        # 如果邻居节点是易感状态，以beta的概率感染邻居节点
                        if G.nodes[neighbor]['state'] == 0 and np.random.random() < beta:
                            G.nodes[neighbor]['state'] = 1
        # 计算当前时间步下的感染节点数，并添加到列表中
        infected_count.append(sum(nx.get_node_attributes(G, 'state').values()))
    # 返回感染节点数列表
    return infected_count

# 定义一些参数和变量
n = 1000 # 节点数
p = 0.01 # 连边概率
beta = 0.2 # 传染率
gamma = 0.1 # 恢复率
T = 100 # 时间步数

# 生成一个随机网络
G = generate_random_network(n, p)

# 模拟SIS模型下的疾病传播
infected_count_SIS = simulate_SIS(G, beta, gamma, T)

# 模拟SIR模型下的疾病传播
infected_count_SIR = simulate_SIR(G, beta, gamma, T)

# 导入matplotlib库
import matplotlib.pyplot as plt

# 绘制SIS模型和SIR模型下的感染节点数随时间变化的曲线图
plt.plot(range(T + 1), infected_count_SIS, label='SIS')
plt.plot(range(T + 1), infected_count_SIR, label='SIR')
plt.xlabel('Time')
plt.ylabel('Infected Count')
plt.legend()
plt.show()