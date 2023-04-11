import numpy as np
import random
import matplotlib.pyplot as plt
import networkx as nx
from pylab import mpl
import csv
mpl.rcParams["font.sans-serif"] = ["SimHei"]
# 定义常量
N = 3000# 网络中节点个数
m = 3 # 每时刻活跃节点发边数量
gama = 2.2 # 定义活跃度幂律分布指数
lBound = 0.001 # 定义活跃度下限
uBound = 1 # 定义活跃度上限                                    /
AVE = 1 # 活跃度+SIS流行病过程平均次数
step = 2000 # 稳态终止时间
eta = 1
infectSeeds = 0.01 # 初始I态节点比例
miuRecover = 0.01 # 定义恢复率，常数,
lanbtaInter = 0.05# 定义
lanbtaInfect0 = 0.055 # 设置初始的感染值
LanbtaMax = 1.5 # 最大感染概率

# 定义活跃度驱动模型类，包含网络构建和流行病传播等方法
class ActivityDrivenModel:
    def __init__(self):
        self.G = nx.Graph() # 创建一个空的无向图对象

    # 根据幂律分布生成活跃度值，并创建节点对象存入图中，同时为每个节点添加活跃度和状态属性
    def generateActive(self):
        for i in range(N):
            activeValue = (uBound ** (1 - gama) - lBound ** (1 - gama)) * random.random() + lBound ** (1 - gama)
            activeValue = activeValue ** (1 / (1 - gama))#
            self.G.add_node(i) # 添加节点到图中
            self.G.nodes[i]['active'] = activeValue # 设置节点的活跃度属性
            self.G.nodes[i]['state'] = 'S' # 设置节点的状态属性

    # 根据活跃度驱动模型构建时序网络，并返回每个时刻的邻接矩阵列表和平均度列表
    def generateRandomNetwork(self):
        adjMatrixList = [] # 存储每个时刻的邻接矩阵的列表
        aveDegreeList = [] # 存储每个时刻的平均度的列表

        for t in range(step):
            # 遍历每个时刻
            self.G.remove_edges_from(self.G.edges())  # 清除图中的所有边
            activeNodeList = []  # 存储每个时刻的活跃节点列表
            for i in range(N):  # 遍历每个节点
                if random.random() < self.G.nodes[i]['active']:  # 如果随机数小于节点的活跃度
                    activeNodeList.append(i)  # 将节点加入到活跃节点列表中
            for i in activeNodeList:  # 遍历每个活跃节点
                neighborList = random.sample(range(N), m)  # 从所有节点中随机选择m个作为邻居节点
                for j in neighborList:  # 遍历每个邻居节点
                    if not self.G.has_edge(i, j):  # 如果图中没有边<i,j>
                        self.G.add_edge(i, j)  # 添加边<i,j>到图中
                        self.G[i][j]['weight'] = 1  # 设置边的权重属性为1
                    else:  # 如果图中已经有边<i,j>
                        self.G[i][j]['weight'] += 1  # 将边的权重属性加1
            adjMatrix = nx.to_numpy_matrix(self.G)  # 将图转换为邻接矩阵
            adjMatrixList.append(adjMatrix)  # 将邻接矩阵加入到列表中
            aveDegree = np.mean(np.sum(adjMatrix, axis=0))  # 计算平均度
            aveDegreeList.append(aveDegree)  # 将平均度加入到列表中

        return adjMatrixList, aveDegreeList

        # 根据SIS模型进行流行病传播，并返回每个时刻的感染比例列表和恢复比例列表
    def SIS(self, adjMatrixList):
        global recoverCount
        infectCountList = []  #每个时刻的感染人数列表
        recoverCountList = [] #每个时刻的恢复人数列表
        infectRatioList = []  # 存储每个时刻的感染比例的列表
        recoverRatioList = []  # 存储每个时刻的恢复比例的列表

        infectSeedList = random.sample(range(N), int(N * infectSeeds))  # 从所有节点中随机选择一定比例作为初始感染种子节点
        for i in infectSeedList:  # 遍历每个感染种子节点
            self.G.nodes[i]['state'] = 'I'  # 将节点状态设置为I

        for t in range(step):  # 遍历每个时刻
            infectCount = 0  # 记录感染节点的数量
            recoverCount = 0  # 记录恢复节点的数量
            for i in range(N):  # 遍历每个节点
                if self.G.nodes[i]['state'] == 'I':  # 如果节点状态是I
                    infectCount += 1  # 感染节点数量加1
                    if random.random() < miuRecover:  # 如果随机数小于恢复率
                        self.G.nodes[i]['state'] = 'S'  # 将节点状态设置为S
                        recoverCount += 1  # 恢复节点数量加1
                elif self.G.nodes[i]['state'] == 'S':  # 如果节点状态是S
                    neighborList = list(self.G.neighbors(i))  # 获取节点的邻居列表
                    for j in neighborList:  # 遍历每个邻居节点
                        if self.G.nodes[j]['state'] == 'I':  # 如果邻居节点状态是I
                            lanbtaInfect = lanbtaInfect0 + lanbtaInter * (self.G[i][j]['weight'] - 1) / (
                                            LanbtaMax - 1) * (
                                                           LanbtaMax - lanbtaInfect0) / LanbtaMax * eta;  # 计算感染率，参考吸引力文章公式（5）
                            if random.random() < lanbtaInfect: # 如果随机数小于感染率
                                self.G.nodes[i]['state'] = 'I' # 将节点状态设置为I
                                break # 跳出循环
            infectCountList.append(infectCount)
            infectRatio = infectCount / N # 计算感染比例
            infectRatioList.append(infectRatio) # 将感染比例加入到列表中

            recoverCountList.append(recoverCount)
            recoverRatio = recoverCount / N # 计算恢复比例
            recoverRatioList.append(recoverRatio) # 将恢复比例加入到列表中

        return infectRatioList, recoverRatioList, infectCountList, recoverCountList

    # 根据SIR模型进行流行病传播，并返回每个时刻的感染比例列表和恢复比例列表
    def SIR(self, adjMatrixList):
        infectCountList = []  # 每个时刻的感染人数列表
        recoverCountList = []  # 每个时刻的恢复人数列表
        infectRatioList = [] # 存储每个时刻的感染比例的列表
        recoverRatioList = [] # 存储每个时刻的恢复比例的列表

        infectSeedList = random.sample(range(N), int(N * infectSeeds)) # 从所有节点中随机选择一定比例作为初始感染种子节点
        for i in infectSeedList: # 遍历每个感染种子节点
            self.G.nodes[i]['state'] = 'I' # 将节点状态设置为I

        for t in range(step): # 遍历每个时刻
            infectCount = 0 # 记录感染节点的数量
            recoverCount = 0 # 记录恢复节点的数量
            for i in range(N): # 遍历每个节点
                if self.G.nodes[i]['state'] == 'I': # 如果节点状态是I
                    infectCount += 1 # 感染节点数量加1
                    if random.random() < miuRecover: # 如果随机数小于恢复率
                        self.G.nodes[i]['state'] = 'R' # 将节点状态设置为R
                        recoverCount += 1 # 恢复节点数量加1
                elif self.G.nodes[i]['state'] == 'S': # 如果节点状态是S
                    neighborList = list(self.G.neighbors(i)) # 获取节点的邻居列表
                    for j in neighborList: # 遍历每个邻居节点
                        if self.G.nodes[j]['state'] == 'I': # 如果邻居节点状态是I
                            lanbtaInfect = lanbtaInfect0 + lanbtaInter * (self.G[i][j]['weight'] - 1) / (LanbtaMax - 1) * (LanbtaMax - lanbtaInfect0) / LanbtaMax * eta;# 计算感染率，参考吸引力文章公式（5）
                            if random.random() < lanbtaInfect: # 如果随机数小于感染率
                                self.G.nodes[i]['state'] = 'I' # 将节点状态设置为I
                                break # 跳出循环
            infectCountList.append(infectCount)
            infectRatio = infectCount / N  # 计算感染比例
            infectRatioList.append(infectRatio)  # 将感染比例加入到列表中

            recoverCountList.append(recoverCount)
            recoverRatio = recoverCount / N  # 计算恢复比例
            recoverRatioList.append(recoverRatio)  # 将恢复比例加入到列表中

        return infectRatioList, recoverRatioList,infectCountList,recoverCountList

    # 根据SEIR模型进行流行病传播，并返回每个时刻的感染比例列表和恢复比例列表
    def SEIR(self, adjMatrixList):
        infectCountList = []  # 每个时刻的感染人数列表
        recoverCountList = []  # 每个时刻的恢复人数列表
        infectRatioList = []  # 存储每个时刻的感染比例的列表
        recoverRatioList = []  # 存储每个时刻的恢复比例的列表

        infectSeedList = random.sample(range(N), int(N * infectSeeds)) # 从所有节点中随机选择一定比例作为初始感染种子节点
        for i in infectSeedList: # 遍历每个感染种子节点
            self.G.nodes[i]['state'] = 'E' # 将节点状态设置为E

        for t in range(step): # 遍历每个时刻
            infectCount = 0 # 记录感染节点的数量
            recoverCount = 0 # 记录恢复节点的数量
            for i in range(N): # 遍历每个节点
                if self.G.nodes[i]['state'] == 'E': # 如果节点状态是E
                    if random.random() < miuRecover: # 如果随机数小于恢复率
                        self.G.nodes[i]['state'] = 'I' # 将节点状态设置为I
                elif self.G.nodes[i]['state'] == 'I': # 如果节点状态是I
                    infectCount += 1 # 感染节点数量加1
                    if random.random() < miuRecover: # 如果随机数小于恢复率
                        self.G.nodes[i]['state'] = 'R' # 将节点状态设置为R
                        recoverCount += 1 # 恢复节点数量加1
                elif self.G.nodes[i]['state'] == 'S': # 如果节点状态是S
                    neighborList = list(self.G.neighbors(i)) # 获取节点的邻居列表
                    for j in neighborList: # 遍历每个邻居节点
                        if self.G.nodes[j]['state'] == 'I': # 如果邻居节点状态是I
                            lanbtaInfect = lanbtaInfect0 + lanbtaInter * (self.G[i][j]['weight'] - 1) / (LanbtaMax - 1) * (LanbtaMax - lanbtaInfect0) / LanbtaMax * eta;# 计算感染率，参考吸引力文章公式（5）
                            if random.random() < lanbtaInfect: # 如果随机数小于感染率
                                self.G.nodes[i]['state'] = 'E' # 将节点状态设置为E
                                break # 跳出循环
            infectCountList.append(infectCount)
            infectRatio = infectCount / N  # 计算感染比例
            infectRatioList.append(infectRatio)  # 将感染比例加入到列表中

            recoverCountList.append(recoverCount)
            recoverRatio = recoverCount / N  # 计算恢复比例
            recoverRatioList.append(recoverRatio)  # 将恢复比例加入到列表中

        return infectRatioList, recoverRatioList,infectCountList,recoverCountList

    # 绘制平均度随时间变化的折线图
    def plotAveDegree(self, aveDegreeList):
        plt.figure() # 创建一个新的图形
        plt.plot(range(step), aveDegreeList) # 绘制折线图
        plt.xlabel("时间") # 设置x轴标签
        plt.ylabel("平均度") # 设置y轴标签
        plt.title("平均度和时间的关系") # 设置标题
        plt.show() # 显示图形

    # 绘制感染比例和恢复比例随时间变化的折线图
    def plotInfectRecover(self, infectRatioList, recoverRatioList):
        plt.figure() # 创建一个新的图形
        plt.plot(range(step), infectRatioList, label="感染率") # 绘制感染比例折线图，并设置标签
        plt.plot(range(step), recoverRatioList, label="恢复率") # 绘制恢复比例折线图，并设置标签
        plt.xlabel("时间") # 设置x轴标签
        plt.ylabel("比例") # 设置y轴标签
        plt.title("感染率和恢复率与时间的关系") # 设置标题
        plt.legend() # 显示图例
        plt.show() # 显示图形
    def plotinfectCount(self,infectCount):
        plt.figure()
        plt.plot(range(step),infectCount,label='感染数')
        plt.xlabel("时间")
        plt.ylabel("感染人数随时间的变化")
        plt.legend()
        plt.show()
    def plotrecoverCount(self,recoverCount):
        plt.figure()
        plt.plot(range(step),recoverCount,label='感染数')
        plt.xlabel("时间")
        plt.ylabel("恢复人数随时间的变化")
        plt.legend()
        plt.show()


# 主函数，测试代码功能
if __name__ == "__main__":
    model = ActivityDrivenModel() # 创建活跃度驱动模型对象
    model.generateActive() # 生成活跃度值并创建节点对象
    adjMatrixList, aveDegreeList = model.generateRandomNetwork() # 构建时序网络并返回邻接矩阵列表和平均度列表
    infectRatioList, recoverRatioList,infectCount,recoverCount = model.SEIR(adjMatrixList) # 根据SIS模型进行流行病传播并返回感染比例列表和恢复比例列表
    model.plotAveDegree(aveDegreeList) # 绘制平均度随时间变化的折线图
    model.plotInfectRecover(infectRatioList, recoverRatioList) # 绘制感染比例和恢复比例随时间变化的折线图
    model.plotinfectCount(infectCount)
    model.plotrecoverCount(recoverCount)