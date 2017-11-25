# This class wraps a variable


class Variable():
    # name is a string
    # values is a list of the possible values for the variable
    def __init__(self, name, values):
        self.name = name
        self.values = values