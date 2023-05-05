import numpy as np
import random
import matplotlib.pyplot as plt
import networkx as nx
from pylab import mpl
import pandas as pd
import csv
mpl.rcParams["font.sans-serif"] = ["SimHei"]
N = 200
m = 7
gama = 2.6
lBound = 0.2
uBound = 5
AVE = 1
step = 500
eta = 1
infectSeeds = 0.1 # 初始I态节点比例
miuRecover = 0.01# 定义恢复率，常数,
lanbtaInter = 0.05 #0.05
lanbtaInfect0 = 1 # 设置初始的感染值
LanbtaMax = 3 # 最大感染概率

# 定义活跃度驱动模型类，包含网络构建和流行病传播等方法
class ActivityDrivenModel:
    def __init__(self):
        self.G = nx.Graph() # 创建一个空的无向图对象

    # 根据幂律分布生成活跃度值，并创建节点对象存入图中，同时为每个节点添加活跃度和状态属性
    def generateActive(self):
        activeList=[]
        for i in range(N):
            activeValue = (uBound ** (1 - gama) - lBound ** (1 - gama)) * random.random() + lBound ** (1 - gama)
            activeValue = activeValue ** (1 / (1 - gama))#
            self.G.add_node(i)
            self.G.nodes[i]['active'] = activeValue
            activeList.append(activeValue)
            self.G.nodes[i]['state'] = 'S'
        return activeList
    # 根据活跃度驱动模型构建时序网络，并返回每个时刻的邻接矩阵列表和平均度列表
    def generateRandomNetwork(self):
        adjMatrixList = [] # 存储每个时刻的邻接矩阵的列表
        aveDegreeList = [] # 存储每个时刻的平均度的列表
        avgPathLengthList=[] #储存每个时刻的平均路径长度
        averageClustList=[]
        for t in range(step):
            # 遍历每个时刻
            self.G.remove_edges_from(self.G.edges())  # 清除图中的所有边
            activeNodeList = []
            for i in range(N):  # 遍历每个节点
                if random.random() < self.G.nodes[i]['active']:
                    activeNodeList.append(i)
            for i in activeNodeList:  # 遍历每个活跃节点
                neighborList = random.sample(range(N), m)
                for j in neighborList:  # 遍历每个邻居节点
                    if not self.G.has_edge(i, j):
                        self.G.add_edge(i, j)
                        self.G[i][j]['weight'] = 1  # 设置边的权重属性为1
                    else:  # 如果图中已经有边<i,j>
                        self.G[i][j]['weight'] += 1
            total_avg_path_length = 0
            total_nodes = 0

            # 遍历G的每个连通分量
            for C in nx.connected_components(self.G):
                # 获取由分量诱导的子图
                S = self.G.subgraph(C)
                # 获取分量中的节点数
                n = S.number_of_nodes()
                # 获取分量中的平均最短路径长度
                l = nx.average_shortest_path_length(S)
                # 按照分量的加权平均值更新总平均路径长度
                total_avg_path_length += n * l
                total_nodes += n

            # 计算最终的平均路径长度，通过除以总节点数
            final_avg_path_length = total_avg_path_length / total_nodes
            avgPathLengthList.append(final_avg_path_length)

            averageClust=nx.average_clustering(self.G)
            averageClustList.append(averageClust)

            adjMatrix = nx.to_numpy_matrix(self.G,dtype='float32')  # 将图转换为邻接矩阵
            adjMatrixList.append(adjMatrix)  # 将邻接矩阵加入到列表中
            aveDegree = np.mean(np.sum(adjMatrix, axis=0))  # 计算平均度
            aveDegreeList.append(aveDegree)  # 将平均度加入到列表中

        return adjMatrixList, aveDegreeList ,avgPathLengthList,averageClustList

        # 根据SIS模型进行流行病传播，并返回每个时刻的感染比例列表和恢复比例列表
    def SI(self, adjMatrixList):
        infectRatioList = []
        infectSeedList = random.sample(range(N), int(N * infectSeeds))  # 从所有节点中随机选择一定比例作为初始感染种子节点
        for i in infectSeedList:
            self.G.nodes[i]['state'] = 'I'
        for t in range(step):
            infectCount = 0
            for i in range(N):
                if self.G.nodes[i]['state'] == 'I':
                    infectCount += 1
                elif self.G.nodes[i]['state'] == 'S':
                    neighborList = list(self.G.neighbors(i))  # 获取节点的邻居列表
                    for j in neighborList:  # 遍历每个邻居节点
                        if self.G.nodes[j]['state'] == 'I':  # 如果邻居节点状态是I
                            lanbtaInfect = lanbtaInfect0 + lanbtaInter * (self.G[i][j]['weight'] - 1) / (
                                            LanbtaMax - 1) * (
                                                           LanbtaMax - lanbtaInfect0) / LanbtaMax * eta;  # 计算感染率，参考吸引力文章公式（5）
                            if random.random() < lanbtaInfect:
                                self.G.nodes[i]['state'] = 'I'
                                break # 跳出循环
            infectRatio = infectCount / N
            infectRatioList.append(infectRatio)

        return infectRatioList

    def SIS(self, adjMatrixList):
        infectRatioList = []
        recoverRatioList = []

        infectSeedList = random.sample(range(N), int(N * infectSeeds))  # 从所有节点中随机选择一定比例作为初始感染种子节点
        for i in infectSeedList:
            self.G.nodes[i]['state'] = 'I'

        for t in range(step):
            infectCount = 0
            recoverCount = 0
            for i in range(N):
                if self.G.nodes[i]['state'] == 'I':
                    infectCount += 1
                    if random.random() < miuRecover:
                        self.G.nodes[i]['state'] = 'S'
                        recoverCount += 1
                elif self.G.nodes[i]['state'] == 'S':
                    neighborList = list(self.G.neighbors(i))  # 获取节点的邻居列表
                    for j in neighborList:  # 遍历每个邻居节点
                        if self.G.nodes[j]['state'] == 'I':  # 如果邻居节点状态是I
                            lanbtaInfect = lanbtaInfect0 + lanbtaInter * (self.G[i][j]['weight'] - 1) / (
                                            LanbtaMax - 1) * (
                                                           LanbtaMax - lanbtaInfect0) / LanbtaMax * eta;  # 计算感染率，参考吸引力文章公式（5）
                            if random.random() < lanbtaInfect:
                                self.G.nodes[i]['state'] = 'I'
                                break # 跳出循环
            infectRatio = infectCount / N
            infectRatioList.append(infectRatio)

            recoverRatio = recoverCount / N
            recoverRatioList.append(recoverRatio)

        return infectRatioList, recoverRatioList

    # 根据SIR模型进行流行病传播，并返回每个时刻的感染比例列表和恢复比例列表
    def SIR(self, adjMatrixList):
        infectRatioList = []
        recoverRatioList = []

        infectSeedList = random.sample(range(N), int(N * infectSeeds))
        for i in infectSeedList:
            self.G.nodes[i]['state'] = 'I'

        for t in range(step):
            infectCount = 0
            recoverCount = 0
            for i in range(N):
                if self.G.nodes[i]['state'] == 'I':
                    infectCount += 1
                    if random.random() < miuRecover:
                        self.G.nodes[i]['state'] = 'R'
                        recoverCount += 1
                elif self.G.nodes[i]['state'] == 'S':
                    neighborList = list(self.G.neighbors(i))
                    for j in neighborList:
                        if self.G.nodes[j]['state'] == 'I':
                            lanbtaInfect = lanbtaInfect0 + lanbtaInter * (self.G[i][j]['weight'] - 1) / (LanbtaMax - 1) * (LanbtaMax - lanbtaInfect0) / LanbtaMax * eta;# 计算感染率，参考吸引力文章公式（5）
                            if random.random() < lanbtaInfect:
                                self.G.nodes[i]['state'] = 'I'
                                break
            infectRatio = infectCount / N
            infectRatioList.append(infectRatio)

            recoverRatio = recoverCount / N
            recoverRatioList.append(recoverRatio)

        return infectRatioList, recoverRatioList
    # 根据SEIR模型进行流行病传播，并返回每个时刻的感染比例列表和恢复比例列表
    def SEIR(self, adjMatrixList):
        infectRatioList = []
        recoverRatioList = []
        exposureRatioList = []
        infectSeedList = random.sample(range(N), int(N * infectSeeds))
        for i in infectSeedList:
            self.G.nodes[i]['state'] = 'E'

        for t in range(step):
            infectCount = 0
            recoverCount = 0
            exposureCount=0
            for i in range(N):
                if self.G.nodes[i]['state'] == 'E':
                    if random.random() < miuRecover:
                        self.G.nodes[i]['state'] = 'I'
                elif self.G.nodes[i]['state'] == 'I':
                    infectCount += 1
                    if random.random() < miuRecover:
                        self.G.nodes[i]['state'] = 'R'
                        recoverCount += 1
                elif self.G.nodes[i]['state'] == 'S':
                    neighborList = list(self.G.neighbors(i))
                    for j in neighborList:
                        if self.G.nodes[j]['state'] == 'I':
                            lanbtaInfect = lanbtaInfect0 + lanbtaInter * (self.G[i][j]['weight'] - 1) / (LanbtaMax - 1) * (LanbtaMax - lanbtaInfect0) / LanbtaMax * eta;# 计算感染率，参考吸引力文章公式（5）
                            if random.random() < lanbtaInfect:
                                self.G.nodes[i]['state'] = 'E'
                                exposureCount +=1
                                break

            infectRatio = infectCount / N
            infectRatioList.append(infectRatio)
            recoverRatio = recoverCount / N
            recoverRatioList.append(recoverRatio)
            exposureRatio=exposureCount/N
            exposureRatioList.append(exposureRatio)


        return infectRatioList, recoverRatioList ,exposureRatioList

# 绘制平均度随时间变化的折线图
    def plotAveDegree(self, aveDegreeList):
        plt.figure()
        plt.plot(range(step), aveDegreeList)
        plt.xlabel("时间")
        plt.ylabel("平均度")
        plt.title("平均度和时间的关系")
        plt.show()
    def plotSI(self,infectRaotioList):
        plt.figure()
        plt.plot(range(step), infectRaotioList, label="感染率")
        plt.xlabel("时间")
        plt.ylabel("比例")
        plt.title("SI感染率和恢复率与时间的关系")
        plt.legend()
        plt.show()
    # 绘制感染比例和恢复比例随时间变化的折线图
    def plotSIR(self, infectRatioList, recoverRatioList):
        plt.figure() # 创建一个新的图形
        plt.plot(range(step), infectRatioList, label="感染率")
        plt.plot(range(step), recoverRatioList, label="恢复率")
        plt.xlabel("时间")
        plt.ylabel("比例")
        plt.title("SIR感染率和恢复率与时间的关系")
        plt.legend()
        plt.show()

    def plotSIS(self, infectRatioList, recoverRatioList):
        plt.figure() # 创建一个新的图形
        plt.plot(range(step), infectRatioList, label="感染率")
        plt.plot(range(step), recoverRatioList, label="恢复率")
        plt.xlabel("时间")
        plt.ylabel("比例")
        plt.title("SIS感染率和恢复率与时间的关系")
        plt.legend()
        plt.show()
    def plotaveragePath(self,averageShortestPathList):
        plt.figure()  # 创建一个新的图形
        plt.plot(range(step), averageShortestPathList, label="平均路径长度")
        plt.xlabel("时间")
        plt.ylabel("平均路径长度")
        plt.title("平均路径长度随时间的变化")
        plt.legend()
        plt.show()

    def plotSIER(self,infectRatioList,recoverRatioList,exposureRatioList):
        plt.figure()  # 创建一个新的图形
        plt.plot(range(step), infectRatioList, label="感染率")
        plt.plot(range(step), recoverRatioList, label="恢复率")
        plt.plot(range(step),exposureRatioList , label="暴露率")
        plt.xlabel("时间")
        plt.ylabel("比例")
        plt.title("SIER感染率和恢复率与时间的关系")
        plt.legend()
        plt.show()
    def plotnxagelustering(self,clusterList):
        plt.figure()  # 创建一个新的图形
        plt.plot(range(step), clusterList, label="平均聚类系数")
        plt.xlabel("时间")
        plt.title("平均聚类系数")
        plt.legend()
        plt.show()

# 主函数，测试代码功能
if __name__ == "__main__":
    model = ActivityDrivenModel() # 创建活跃度驱动模型对象
    activelist=model.generateActive() # 生成活跃度值并创建节点对象
    data2 = pd.DataFrame(data=activelist, index=None)
    data2.to_csv('activeValue')
    adjMatrixList, aveDegreeList ,avgPathLengthList,agelustering= model.generateRandomNetwork() # 构建时序网络并返回邻接矩阵列表和平均度列表
    model.plotAveDegree(aveDegreeList) # 绘制平均度随时间变化的折线图
    model.plotaveragePath(avgPathLengthList)
    model.plotnxagelustering(agelustering)

    # infectRatioList=model.SI(adjMatrixList)
    # model.plotSI(infectRatioList)

    # infectRatioList, recoverRatioList = model.SIS(adjMatrixList) # 根据SIS模型进行流行病传播并返回感染比例列表和恢复比例列表
    # model.plotSIS(infectRatioList, recoverRatioList) # 绘制感染比例和恢复比例随时间变化的折线图
    #
    # infectRatioList3, recoverRatioList3 = model.SIR(adjMatrixList)  # 根据SIS模型进行流行病传播并返回感染比例列表和恢复比例列表
    # model.plotSIR(infectRatioList3, recoverRatioList3)  # 绘制感染比例和恢复比例随时间变化的折线图
    #
    infectRatioList4, recoverRatioList4, exposureRatioList4 = model.SEIR(adjMatrixList)  # 根据SIS模型进行流行病传播并返回感染比例列表和恢复比例列表
    model.plotSIER(infectRatioList4, recoverRatioList4,exposureRatioList4)  # 绘制感染比例和恢复比例随时间变化的折线图



