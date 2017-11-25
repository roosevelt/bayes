# A conditional probability distribution table to be associated to a node
from decimal import Decimal
import copy

class CPD():
    # variables is a list of variables of type Variable
    # parameters is a list of 
    def __init__(self, variable, condition_variables=[]):
        if len(condition_variables)>0:
            condition_variables.sort(key = lambda x: x.name)
        self.variable = variable
        self.condition_variables = condition_variables
        self.table = self.__buildTable(self.variable, self.condition_variables) 
        
    def __buildTable(self, variable, condition_variables):
        table = {}
        for combination in self.__getVariableValuesCombinations([variable]+condition_variables):
            f_comb = frozenset(combination)
            table[f_comb] = Decimal('0.0')
        return table
    
    def __getVariableValuesCombinations(self, variables):
        # Generate all possible combinations of data
        variable_values = {}
        variable_card_prod = {}
        count = 1
        
        for variable in variables:
            variable_values[variable.name] = variable.values
            variable_card_prod[variable.name] = count
            count = count * len(variable.values)
        
        #Generate all variables combinations probabilities
        var_counters = {}
        var_counters = copy.copy(variable_card_prod)
        var_actual_values = [0] * len(variables)
        k = 0
        while (k < count):
            example = []
            for i in range(len(variables)):
                var_name = variables[i].name
                value = variable_values[var_name][var_actual_values[i]]
                var_counters[var_name] = var_counters[var_name] - 1
                if var_counters[var_name] == 0:
                    var_counters[var_name] = variable_card_prod[var_name]
                    if (var_actual_values[i] != ((len(variable_values[var_name])-1))):
                        var_actual_values[i] = var_actual_values[i] + 1
                    else:
                        var_actual_values[i] = 0
                example.append((var_name,value))
            yield example
            k+=1;
        
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
        result += "Lines: " +  str(len(self.table)) +" entries\n\n"
    
        for example in self.table:
           str_example = ""
           feature_val_lst = []
           main_feature = ""
           for feature_val in example:
               if(feature_val[0]!=self.variable.name):
                   feature_val_lst.append(feature_val[0] + "(" + str(feature_val[1]) + ")")
                   feature_val_lst.sort()
               else:
                   main_feature = feature_val[0] + "(" + str(feature_val[1]) + ")"
           feature_val_lst = [main_feature] + feature_val_lst

           for feat_str in feature_val_lst:
               str_example += '{0:30s} '.format(feat_str)
           str_example += '{0:25s} '.format(str(self.table[example]))
           result+=str_example+"\n"
        return result
    
    # Set the parameters in the CPD
    def setParameters(self, parameters):
        # parameters is a list of parameters in the order that it appears on the table
        k=0
        for combination in self.table:
            self.table[combination]=Decimal(str(parameters[k]))
            k+=1
        
    def getParameters(self):
        return self.table