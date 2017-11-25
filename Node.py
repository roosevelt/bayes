# A Bayesian Node, with its properties and associated CPD
class Node():
    
    def __init__(self, CPD):
        self.name = CPD.variable.name
        self.CPD = CPD
        self.parents = CPD.condition_variables
        
    def __str__(self):
        return self.CPD.__str__()
        
    
        