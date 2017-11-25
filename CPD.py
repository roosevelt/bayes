# A conditional probability distribution table to be associated to a node
import decimal
import copy

class CPD():
    # variables is a list of variables of type Variable
    # parameters is a list of 
    def __init__(self, variables):
        variables.sort(key = lambda x: x.name)
        self.variables = variables
        self.table = self.__buildTable(self.variables) 
    
    def setParameters(parameters):
        pass
    
    def __buildTable(self, variables):
        table = {}
        for combination in self.__getVariableValuesCombinations(variables):
            f_comb = frozenset(combination)
            table[f_comb] = decimal.Decimal('0.0')
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
        result = ""
        for example in self.table:
           str_example = ""
           feature_val_lst = []

           for feature_val in example:
               feature_val_lst.append(feature_val[0] + "(" + str(feature_val[1]) + ")")
               feature_val_lst.sort()

           for feat_str in feature_val_lst:
               str_example += '{0:30s} '.format(feat_str)
           str_example += '{0:25s} '.format(str(self.table[example]))
           result+=str_example+"\n"
        return result