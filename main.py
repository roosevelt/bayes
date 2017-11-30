from Variable import Variable
from CPD import CPD
from Node import Node
from BN import BayesianNetwork
from Inference import Inference
 
bn2 = BayesianNetwork()
bn2.loadBIF('./data/asia.bif')
print(bn2)

#inf = Inference()
#print([(node.name, node.variable.values) for node in bn2.nodes])
#print(inf.variableElimination([('DISCONNECT', 'FALSE'),('ERRLOWOUTPUT', 'FALSE'), ('LVFAILURE', 'FALSE'), ('KINKEDTUBE', 'FALSE')],[], bn2))



#inference_engine = Inference()
#result = inference_engine.variableElimination(bn2, query=[('asia', 'yes')])
#print(result)

