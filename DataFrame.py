import numpy as np

class MCDataFrame:
    def __init__(self, filename):
        with open(filename, 'r') as f:
            lines = f.readlines()

            first_line = lines[0]
            self.data = {}
            for n,line in enumerate(lines[1::2]):
                self.data[lines[n+1]] = lines[n+2]
                


d = MCDataFrame('testing.txt')
