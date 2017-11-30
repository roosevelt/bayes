# A conditional probability distribution table to be associated to a node
from Table import Table

class CPD():
    # variables is a list of variables of type Variable
    # parameters is a list of 
    def __init__(self, variable, condition_variables=[]):
        if len(condition_variables)>0:
            condition_variables.sort(key = lambda x: x.name)
        self.variable = variable
        self.condition_variables = condition_variables
        self.table = Table([self.variable]+self.condition_variables) 
        
        self.__variable_dict = self.__buildVariableDict(self.variable, self.condition_variables)
        
    def __str__(self):
        # Prints this CPD to the stdout
        n_line = 20
        result = ""
        result += ("="*n_line) + " Conditional Probability Table of " + self.variable.name  + " (" + str(len(self.variable.values))+ " Values) " + "="*n_line + "\n"
        condition_arr = []
        for condition_var in self.condition_variables:
            condition_arr.append(condition_var.name + " ("  + str(len(condition_var.values)) +  " Values)")        
        result += "Possible Values: " + ", ".join(self.variable.values) + "\n"
        if len(self.condition_variables)>0:
            result +=  "Parents: " + ", ".join(condition_arr) + "\n"
        else:
            result += "Parents: None (no parents)\n"
        result += "Lines: " +  str(len(self.table.data)) +" entries\n\n"

        self.table.setFirstVariableToPrint(self.variable)
        result += self.table.__str__()
        
        return result
    
    # Set the parameters in the CPD
    def setParameters(self, parameters):
        # parameters is a list of parameters in the order that it appears on the table
        self.table.setParameters(parameters)
        
    def getParameters(self):
        return self.table.getParameters()
    
    def getVariableByName(self, name):
        return self.__variable_dict[name]
        
    def __buildVariableDict(self, variable, conditional_variable):
        variables = [variable] + conditional_variable
        variable_dict = {}
        for variable in variables:
            variable_dict[variable.name] = variable
        return variable_dict