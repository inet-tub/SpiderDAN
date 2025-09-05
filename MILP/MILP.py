class Variables:
  def __init__(self,inst,curNumNode):
    self.model = gp.Model("Matching" + str(curNumNode))
    self.active = {} #a_{u,v}
    self.dis = {} #dis_{u,v}
    self.valid = {} # y^w_{u,v}
    for firstNode in inst.nodes:
      self.dis[firstNode] = {}
      self.valid[firstNode] = {}
      self.active[firstNode] = {}
      for secondNode in inst.nodes:
        self.valid[firstNode][secondNode] = {}
        self.dis[firstNode][secondNode] = self.model.addVar(vtype = GRB.INTEGER, name = "dist " +str(firstNode)+"  "+str(secondNode))
        self.active[firstNode][secondNode] = self.model.addVar(vtype = GRB.BINARY, name = "active "+str(firstNode)+"  "+str(secondNode))
        for thirdNode in inst.nodes:
          self.valid[firstNode][secondNode][thirdNode] = self.model.addVar(vtype = GRB.BINARY)
  def construct_constraints(self,inst):
    for edge in inst.edges:
      u,v=edge
      self.model.addLConstr(self.dis[u][v], GRB.EQUAL, 1)
      self.model.addLConstr(self.dis[v][u], GRB.EQUAL, 1)
    for u in inst.nodes:
      SumEdges = gp.LinExpr()
      for v in inst.nodes:
        if (u !=v and inst.notAnEdge(u,v)):
          self.model.addLConstr(self.active[u][v], GRB.EQUAL, self.active[v][u])
          self.model.addLConstr(self.dis[u][v], GRB.GREATER_EQUAL, 1)
          DisActive = gp.LinExpr()
          DisActive.addTerms(1,self.active[u][v])
          DisActive.addTerms(-inst.MaxValue,self.active[u][v])
          self.model.addLConstr(self.dis[u][v], GRB.LESS_EQUAL, DisActive+inst.MaxValue)
          ActivitySum = gp.LinExpr()
          for w in inst.nodes:
            if (w != u and w !=v):
              DisSum = gp.LinExpr()
              DisSum.addTerms(1,self.dis[u][w])
              DisSum.addTerms(1,self.dis[w][v])
              self.model.addLConstr(self.dis[u][v], GRB.LESS_EQUAL, DisSum)
              DisMax = DisSum
              DisMax.addTerms(inst.MaxValue,self.valid[u][v][w])
              self.model.addLConstr(self.dis[u][v], GRB.GREATER_EQUAL, DisMax-inst.MaxValue)
              ActivitySum.addTerms(1,self.valid[u][v][w])
          ActivitySum.addTerms(1,self.active[u][v])
          self.model.addLConstr(ActivitySum, GRB.EQUAL, 1 )
          SumEdges.addTerms(1,self.active[u][v])
      self.model.addLConstr(SumEdges,GRB.EQUAL,numMatching)
      self.model.addLConstr(self.dis[u][u], GRB.EQUAL, 0)
  def construct_objectives(self,inst):
      obj_expr = gp.LinExpr()
      for u in inst.nodes:
        for v in inst.nodes:
          obj_expr.add(int(inst.demand[u][v])*self.dis[u][v])
      self.model.setObjective(obj_expr, GRB.MINIMIZE)

"""## MIP Solver

"""

def MIPSolve(inst,curNumNode):
  vars = Variables(inst,curNumNode)
  vars.construct_constraints(inst)
  vars.construct_objectives(inst)
  vars.model.setParam('OutputFlag', False)   #Quite Opetimzation
  vars.model.optimize()

  if (vars.model.SolCount > 0):
    allVars = vars.model.getVars()
    MIPGraph = nx.Graph()
    return [vars.model.ObjVal, MIPGraph]
  else:
    return -1