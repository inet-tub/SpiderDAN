def findCost(graph,demand):
  allLength = dict(nx.all_pairs_shortest_path_length(graph))
  sumCost = 0
  for v in graph.nodes:
    for u in graph.nodes:
      sumCost+= allLength[u][v]*demand[u][v]
  return sumCost

"""## SuperNode Heuristics

### SuperNode Transoformations
"""

def getDemandDegree(inst):
  count = 0
  average = sum(sum(row) for row in inst.demand) / (len(inst.demand) * len(inst.demand[0]))
  for u in inst.nodes:
    for v in inst.nodes:
      if (inst.demand[u][v] > average):
        count+=1
  return float(count/inst.numNodes)

def superFather(myNode, fact):
  return int(myNode/fact)

def createSuperDemand(numNodes, oldInst, fact): #TODO Arash: Merge demands of consecutive factor nodes
  matrix = [[0] * (numNodes+1) for _ in range(numNodes+1)]
  for u in oldInst.nodes:
    for v in oldInst.nodes:
      matrix[superFather(u,fact)][superFather(v,fact)]+= oldInst.demand[u][v]
  return matrix

def createSuperInstance(inst,demDeg):
  superFactor = math.floor(demDeg)
  superNumNode = math.floor(inst.numNodes/superFactor)
  superDemand = createSuperDemand(superNumNode, inst, superFactor)
  return superDemand, superNumNode

def createSuperInstanceChord(inst,demDeg):
  superFactor = math.floor(demDeg)
  superNumNode = math.floor(inst.numNodes/superFactor)
  return superNumNode

def transferSupertoNormal(oldGraph,matching,factor, nowDemand):
  currGraph = oldGraph.copy()
  flags = [False] * len(oldGraph.nodes)
  for edge in matching:
    u , v = edge
    uChild = getNextChild(flags,u,factor)
    vChild = getNextChild(flags,v,factor)
    if(not((uChild,vChild) in currGraph.edges) and not((vChild,uChild) in currGraph.edges)):
      flags[uChild] = True
      flags[vChild] = True
      currGraph.add_edge(uChild,vChild)
  end_time = time.time()

  demandGraph = buildDemandGraphMatching(nowDemand,len(oldGraph.nodes), oldGraph.edges)
  for edge in currGraph.edges:
    if(edge in demandGraph.edges):
      demandGraph.remove_edge(edge[0],edge[1])
  for v in currGraph.nodes:
    if (flags[v] == True):
      demandGraph.remove_node(v)

  nowMatching = nx.max_weight_matching(demandGraph, maxcardinality=True)
  for edge in nowMatching:
    u,v = edge
    currGraph.add_edge(u,v)
    flags[u] = True
    flags[v] = True

  return currGraph,  end_time

def getNextChildForSuper(curr, fact, numNodes):
  startSimpleNode = curr*fact
  allCur = []
  for i in range(fact):
    currNode = (startSimpleNode + i)% numNodes
    allCur.append(currNode)
  return allCur

def getNextChild(flags, curr, fact):
  startSimpleNode = curr*fact
  for i in range(fact):
    currNode = (startSimpleNode + i)% len(flags)
    if (flags[currNode] == False):
      return currNode
  return "ERROR"


def SpiderTransferSupertoNormal(oldGraph,matching,factor, nowDemand):
  currGraph = oldGraph.copy()
  flags = [False] * len(oldGraph.nodes)
  tempGraph = nx.Graph()
  for node in oldGraph.nodes:
    tempGraph.add_node(node)
  for edge in matching:
    u , v = edge
    uChilds = getNextChildForSuper(u,factor,len(oldGraph))
    vChilds = getNextChildForSuper(v,factor,len(oldGraph))
    for uChild in uChilds:
      for vChild in vChilds:
        tempGraph.add_edge(uChild,vChild, weight = nowDemand[uChild][vChild])
    # if(not((uChild,vChild) in currGraph.edges) and not((vChild,uChild) in currGraph.edges)):
    #   flags[uChild] = True
    #   flags[vChild] = True
    #   currGraph.add_edge(uChild,vChild)
  tempMatching = nx.max_weight_matching(tempGraph, maxcardinality=True)
  for edge in tempMatching:
    u , v = edge
    currGraph.add_edge(u,v)

  end_time = time.time()


  demandGraph = buildDemandGraphMatching(nowDemand,len(oldGraph.nodes), oldGraph.edges)
  for edge in currGraph.edges:
    if(edge in demandGraph.edges):
      demandGraph.remove_edge(edge[0],edge[1])
  for v in currGraph.nodes:
    if (flags[v] == True):
      demandGraph.remove_node(v)

  nowMatching = nx.max_weight_matching(demandGraph, maxcardinality=True)
  for edge in nowMatching:
    u,v = edge
    currGraph.add_edge(u,v)
    flags[u] = True
    flags[v] = True


  return currGraph,  end_time

"""### SuperNode Chord

TODO:

*   Use lambda function instead
*   Fix non power two edges


"""

from scipy.special import lambertw

def createChord(numNodes): #TODO Arash
  myLogN = int(math.log2(numNodes))
  chorded = nx.Graph()
  numNodes= numNodes-1
  for v in range(numNodes):
    for i in range(myLogN):
      nextNode = int((v + math.pow(2,i)) % numNodes)
      chorded.add_edge(v,nextNode)
  return chorded

def prevPowerTwo(curN):
  prevPowTwo = 2 ** int(np.log2(curN))
  return prevPowTwo

def superChord(inst):
  start_time = time.time()
  #inst = transformToPowerTwo(inst)
  #logN = 2*int(lambertw(inst.numNodes*math.log(2))/math.log(2))
  logN = 2*int(math.log2(inst.numNodes))-1
  #print(logN)
  superNumNode = createSuperInstanceChord(inst,logN)
  superNumNode = prevPowerTwo(superNumNode)
  chordMatching = createChord(superNumNode)
  graphAndMatching, end_time  = transferSupertoNormal(inst.graph,chordMatching.edges,logN, inst.demand)
  return findCost(graphAndMatching, inst.demand)

"""### SuperNode ConstantDegree

#### Aleksander's previous code
"""

def get_num_missing_huffman(n, arity=2, fat_root=False):
    while True:
        if fat_root and n <= arity+1:
            return arity+1-n
        if not fat_root and n <= arity:
            return arity-n
        n -= arity
        n += 1


# create huffman tree for D[v] and place it in G
# also save a mapping from edge to leaf in L
def huffman_tree(D, v, G, L, arity=2, fat_root=False):
    #assert(v not in G or G.degree(v) == 0)

    Q = PriorityQueue()

    if v not in L:
        L[v] = dict()
    for u in D.neighbors(v):
        w = D[v][u]["weight"]
        leaf = "T_"+str(v)+"_"+str(u)
        G.add_node(leaf)
        L[v][u] = leaf
        Q.put((w, leaf))

    missing = get_num_missing_huffman(D.degree(v), arity=arity, fat_root=fat_root)
    missing_V = []
    for i in range(missing):
        leaf = "missing_"+str(i)
        G.add_node(leaf)
        Q.put((0, leaf))
        missing_V.append(leaf)

    root = None
    while Q.qsize() > 1:
        wsum = 0
        parent = "huffman_"+str(len(G))
        G.add_node(parent)


        assert(Q.qsize() >= arity)


        merged = []
        for i in range(arity):
            w, node = Q.get()
            wsum += w
            G.nodes[node]["parent"] = parent
            G.add_edge(node, parent)
            merged.append(node)
        if fat_root and Q.qsize() == 1:
            w, node = Q.get()
            wsum += w
            G.nodes[node]["parent"] = parent
            G.add_edge(node, parent)
            merged.append(node)
        if not fat_root:
            assert(G.degree(parent) <= arity)
        else:
            assert(G.degree(parent) <= arity+1)

        if Q.qsize() == 0:
            root = parent
        Q.put((wsum, parent))

    if root is None:
        assert(D.degree(v) == 1)
        root = "huffman_"+str(len(G))
        G.add_node(root)
        u = L[v][list(D.neighbors(v))[0]]
        G.nodes[u]["parent"] = root
        G.add_edge(u, root)

    for m in missing_V:
        G.remove_node(m)

    # relabel root to v
    root_neighbours = list(G.neighbors(root))
    G.remove_node(root)
    G.add_node(v)
    G.nodes[v]["parent"] = None
    for neigh in root_neighbours:
        G.nodes[neigh]["parent"] = v
        G.add_edge(neigh, v)
    L[v][v] = v

    if not fat_root:
        assert(G.degree(v) <= arity)
    else:
        assert(G.degree(v) <= arity+1)

def additive_model(D, max_deg, skip_deg0=False, fat_root=False):
    G = nx.Graph()
    L = dict()
    if not skip_deg0:
        for v in D.nodes:
            G.add_node(v)
    for v in D.nodes:
        if skip_deg0 and D.degree(v) == 0:
            continue
        huffman_tree(D, v, G, L, arity=max_deg-1, fat_root=fat_root)
    for u, v in D.edges():
        Tuv = L[u][v]
        Tvu = L[v][u]

        G.add_edge(G.nodes[Tuv]["parent"], G.nodes[Tvu]["parent"])
        G.remove_node(Tuv)
        G.remove_node(Tvu)

    return G

def rdrg_is_suitable(G, potential_edges):
    for a in potential_edges:
        for b in potential_edges:
            if a == b:
                continue
            if not G.has_edge(a, b):
                return True
    return False

def rdrg_get_stubs(potential_edges):
    stubs = []
    for v, count in potential_edges.items():
        for _ in range(count):
            stubs.append(v)
    random.shuffle(stubs)
    return stubs

def rdrg(G, max_deg):
    potential_edges = dict()
    for v in G.nodes:
        missing = max_deg - G.degree(v)
        if missing > 0:
            potential_edges[v] = missing

    stubs = rdrg_get_stubs(potential_edges)
    while len(stubs) > 1:
        potential_edges = dict()

        for i in range(0, len(stubs), 2):
            j = i + 1
            if j == len(stubs):
                break
            a = stubs[i]
            b = stubs[j]

            if a != b and not G.has_edge(a, b):
                G.add_edge(a, b)
            else:
                for x in [a,b]:
                    if x not in potential_edges:
                        potential_edges[x] = 0
                    potential_edges[x] += 1

        if not rdrg_is_suitable(G, potential_edges):
            break
        stubs = rdrg_get_stubs(potential_edges)


def fixed_degree_model(D, max_deg):
    G = None
    assert(max_deg >= 6)
    d1 = max_deg - 3
    d2 = 3

    wE = []
    for u,v,w in D.edges.data("weight"):
        wE.append((w,u,v))
    wE = list(sorted(wE, reverse=True))


    e_lb = 0
    e_ub = len(wE)

    while e_lb != e_ub:
        e_target = (e_lb + e_ub + 1) // 2
        H = nx.Graph()
        H.add_nodes_from(D)
        for w, u, v in wE[:e_target]:
            H.add_edge(u,v, weight=w)
        G = additive_model(H, d1, skip_deg0=True, fat_root=True)
        success = len(G.nodes) <= len(D.nodes)
        if success:
            e_lb = e_target
        else:
            e_ub = e_target-1

    H = nx.Graph()
    H.add_nodes_from(D)
    G = additive_model(H, d1, skip_deg0=True, fat_root=True)
    for w, u, v in wE[:e_ub]:
        H.add_edge(u,v, weight=w)
    G = additive_model(H, d1, skip_deg0=True, fat_root=True)
    assert(len(G.nodes) <= len(D.nodes))

    # fill remaining degrees
    missing = len(D.nodes) - len(G.nodes) # missing nodes
    for i in range(missing):
        G.add_node("fixed_deg_missing_"+str(i))
    rdrg(G, max_deg)

#    # relabel graph (map huffman internal nodes to nodes in D)
    New = set()
    Vp = set()
    for v in G.nodes():
        if v in D.nodes():
            Vp.add(v)
        else:
            New.add(v)
    Vpp = set(D.nodes) - Vp
    M = dict()
    assert(len(Vpp) == len(New))

    Vpp = list(Vpp)
    New = list(New)
    Vp  = list(Vp)

    for i, v in enumerate(New):
        M[v] = Vpp[i]
    for v in Vp:
        M[v] = v

    G = nx.relabel_nodes(G, M)

    return G

def thresh_balance(D):
    G = nx.Graph()

    n = len(D.nodes())
    m = len(D.edges())

    # average degree, rounded down
    avg_deg = (2 * m) // n
    # 2*average degree, rounded up
    hdeg = (4 * m + n-1) // n
    if hdeg < 2:
        hdeg = 2


    # split graph into high and low degree vertices
    H = set()
    L = set()
    for v in D.nodes():
        if D.degree(v) > hdeg:
            H.add(v)
        else:
            L.add(v)

    # build huffman trees for high deg vertices
    internal = set()

    Lmap = dict()
    for v in H:

        # find low demand edges
        edges = []
        for u in D.neighbors(v):
            edges.append((D[v][u]["weight"], u))
        low_edges = sorted(edges)[:-hdeg]

        D_tmp = nx.Graph()
        for w, u in low_edges:
            D_tmp.add_edge(u, v, weight=w)

        # create huffman tree for low demand edges
        huffman_tree(D_tmp, v, G, Lmap, arity=hdeg)
        assert(G.degree(v) <= hdeg)

        # save internal vertices
        visited = set()
        visited.add(v)

        layer = [v]
        nextlayer = []

        while True:
            for a in layer:
                for b in G.neighbors(a):
                    if b not in visited:
                        if G.degree(b) != 1:
                            internal.add(b)
                        visited.add(b)
                        nextlayer.append(b)
            layer = nextlayer
            nextlayer = []

            if len(layer) == 0:
                break

    for v in G.nodes():
        assert(G.degree(v) <= hdeg + 1)

    # map edges of D

    high_weight = dict()
    for v in D.nodes():
        high_weight[v] = 0

    for u, v in D.edges():
        # is this a low demand edge for u?
        low_u = False if u in L else v in Lmap[u]
        # is this a low demand edge for u?
        low_v = False if v in L else u in Lmap[v]

        a = None
        b = None
        if not low_u:
            high_weight[u] += 1
            a = u
        else:
            a = list(G.neighbors(Lmap[u][v]))[0] # get the parent
            #assert(a in internal)
            G.remove_node(Lmap[u][v])
        if not low_v:
            high_weight[v] += 1
            b = v
        else:
            b = list(G.neighbors(Lmap[v][u]))[0] # get the parent
            #assert(b in internal)
            G.remove_node(Lmap[v][u])
        G.add_edge(a,b)

    for v in G.nodes():
        assert(G.degree(v) <= 2*hdeg + 1)

    # map internal nodes to L nodes
    Imap = dict()

    internal_arr = list(internal)
    L_arr = list(L)
    assert(len(internal_arr) <= len(L))

    for i in range(len(internal_arr)):
        a = internal_arr[i]
        b = L_arr[i]

        Imap[a]=b
    G = nx.relabel_nodes(G, Imap)
    return G

"""#### My code"""

def buildDemandGraph(demand, numNodes, edgeOfGraph):
  demandGraph = nx.Graph()
  for i in range(numNodes):
    for j in range(i+1, numNodes):
      if(demand[i][j] > 0 and not((i,j) in edgeOfGraph)):
        demandGraph.add_edge(i,j, weight = demand[i][j])
  return demandGraph

def buildDemandGraphMatching(demand, numNodes, edgeOfGraph):
  demandGraph = nx.Graph()
  for i in range(numNodes):
    for j in range(numNodes):
      # x = demand[i][j]
      # if (x == 0.0):
      if(not((i,j) in edgeOfGraph)):
        demandGraph.add_edge(i,j, weight = demand[i][j])
  return demandGraph

def findConstantNet(superDemand, superNumNode, degreeOnTop, infEdges):
  demandGraph = buildDemandGraph(superDemand, superNumNode, infEdges)
  return fixed_degree_model(demandGraph, degreeOnTop)

#def transferSupertoNormal(oldGraph,matching,factor):
def superPlusConstDegree(inst):
  #print(inst.numNodes)
  start_time = time.time()
  demDeg = max(6,math.ceil(getDemandDegree(inst)))
  #print(demDeg)
  superDemand, superNumNode = createSuperInstance(inst,demDeg)
  #print(demDeg, superNumNode)
  superMatching = findConstantNet(superDemand, superNumNode, demDeg, inst.edges)
  #nx.draw(superMatching)
  #nx.draw(inst.graph)
  matchings, end_time = SpiderTransferSupertoNormal(inst.graph,superMatching.edges,int(demDeg), inst.demand)
  #nx.draw(matchings)
  return findCost(matchings,inst.demand)

"""## Other Greedy Heuristics

### Maximum matchings
"""

def maxMatchWithDemand(inst):
  demandGraph = buildDemandGraph(inst.demand,inst.numNodes, inst.edges)
  nowMatching = nx.max_weight_matching(demandGraph)
  curGraph = inst.graph.copy()
  for edge in nowMatching:
    u,v = edge
    curGraph.add_edge(u,v)
  return findCost(curGraph,inst.demand)

def maxMatchWithCostReduction(inst):
  costGraph = nx.Graph()
  flags = [False] * inst.numNodes
  curNodes = list(inst.nodes).copy()
  curGraph = inst.graph.copy()
  allCost = findCost(curGraph,inst.demand)
  for u in curNodes:
      for v in curNodes:
        nowGraph = curGraph.copy()
        nowGraph.add_edge(u,v)
        tempCost = findCost(nowGraph,inst.demand)
        costGraph.add_edge(u,v,weight = allCost-tempCost)
  nowMatching = nx.max_weight_matching(costGraph)
  for edge in nowMatching:
    u,v = edge
    curGraph.add_edge(u,v)
  return findCost(curGraph,inst.demand)

"""### Slow greedy"""

#ToDo add flags to ensure that you do not repeat one thing more than once

def addOneReduceTheCostMost(inst):
  flags = [False] * inst.numNodes
  curNodes = list(inst.nodes).copy()
  curGraph = inst.graph.copy()
  for i in range(int(inst.numNodes/2)):
    lowestCost = findCost(curGraph,inst.demand)
    reslutGraph = curGraph.copy()
    tempV = 0
    tempU = 0
    for u in curNodes:
      if (flags[u] == False):
        for v in curNodes:
          if (flags[v] == False):
            nowGraph = curGraph.copy()
            nowGraph.add_edge(u,v)
            tempCost = findCost(nowGraph,inst.demand)
            if ( tempCost < lowestCost):
              resultGraph = nowGraph
              lowestCost = tempCost
              tempV = v
              tempU = u
    flags[tempV] = True
    flags[tempU] = True
    curGraph = resultGraph
  return lowestCost
#addOneReduceTheCostMost(testInst)

"""### Fast greedy"""

def fastGreedy(inst):
  sortedDemand = [(val, i, j) for i, row in enumerate(inst.demand) for j, val in enumerate(row)]
  sortedDemand.sort(reverse=True)
  curGraph = inst.graph.copy()
  flags = [False] * inst.numNodes
  for edge in sortedDemand:
    val, u, v = edge
    if (flags[u]==False and flags[v]==False and u!=v and not((u,v)in curGraph.edges) and not((v,u)in curGraph.edges) ):
      curGraph.add_edge(u,v)
      flags[u] = True
      flags[v] = True
  return findCost(curGraph,inst.demand)

"""

# Running Algs"""

def testAll(infType, demandType, demand_edges, numTests, graphNumNodes):
  rows = 6
  result = [0] * rows
  approx = [0] * rows
  timeAlg = [0] * rows
  for aTest in range(numTests):
    testInst = Instance(infType, demandType, demand_edges, graphNumNodes)
    sizeGraph = testInst.numNodes

    start_time = time.time()
    result[5]+= findCost(testInst.getGraph(),testInst.getDemand())
    end_time = time.time()
    timeAlg[5]+= end_time-start_time

    # print("fastGreedy")

    start_time = time.time()
    result[1]+= fastGreedy(testInst)
    end_time = time.time()
    timeAlg[1]+= end_time-start_time


    # print("maxMatchWithDemand")

    start_time = time.time()
    result[2]+= maxMatchWithDemand(testInst)
    end_time = time.time()
    timeAlg[2]+=end_time-start_time


    # print("superChord")

    start_time = time.time()
    result[3] += superChord(testInst)
    end_time = time.time()
    timeAlg[3] +=end_time-start_time

    # print("superPlusConstDegree")

    start_time = time.time()
    result[4] += superPlusConstDegree(testInst)
    end_time = time.time()
    timeAlg[4] +=end_time-start_time

    # start_time = time.time()
    # MIPRes, MIPGraph = MIPSolve(testInst,sizeGraph)
    # result[0]+=MIPRes
    # end_time = time.time()
    # timeAlg[0]+=end_time-start_time

  for row in range(rows):
      result[row]/=numTests*sizeGraph*(sizeGraph-1)/2 ####### TODO: Normalize results ######
  for row in range(rows):
      timeAlg[row]/=numTests
  return result, timeAlg

"""## Approximation and Absoulte Plots"""

def plotVariables(result, Variables, nameOfVar):

  plot_list = list(map(list, zip(*result)))
  line_names = [""]*timeSize
  line_names[0] = "Fast Greedy"
  line_names[1] = "Slow Greedy"
  line_names[2] = "Matching on demands"
  line_names[3] = "Matching on cost reduction"
  line_names[4] = "Supernodes + Chord"
  #line_names[5] = "MIP"

  for i in range(timeSize):
      plt.plot(plot_list[i], label=line_names[i],marker=Markers[i], color=Colors[i])

  plt.xticks(np.arange(len(Variables)), Variables)
  plt.xlabel(nameOfVar)
  plt.ylabel("Average weighted distance")
  plt.legend()

  plt.show()