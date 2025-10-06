import numpy as np
from binarytree import build
from sys import exit

# ------------------------------------------------------------
class Region:
    '''A binary-splittable 2D-region.'''

    def __init__(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def __repr__(self):
        return f'Region({self.xmin}, {self.xmax}, {self.ymin}, {self.ymax})'

    @property
    def limits(self):
        return self.xmin, self.xmax, self.ymin, self.ymax

    def xsplit(self, value):
        '''Split this region in the x-axis.'''
        r1 = Region(self.xmin, value, self.ymin, self.ymax)
        r2 = Region(value, self.xmax, self.ymin, self.ymax)
        return r1, r2

    def ysplit(self, value):
        '''Split this region in the y-axis.'''
        r1 = Region(self.xmin, self.xmax, self.ymin, value)
        r2 = Region(self.xmin, self.xmax, value, self.ymax)
        return r1, r2

    def split(self, value, axis):
        '''Split this region along the given axis.'''
        if axis == 0:
            return self.xsplit(value)
        elif axis == 1:
            return self.ysplit(value)


# ------------------------------------------------------------
class FindNodes:
    def __init__(self, data):
        self.data = data
        self.nodes = []
        self.n_in_cells = []
        self.book = {}
        self.col = False
        self.level = 1

    def __repr__(self):
        return self.nodes

    def sort_data(self, data, col):
        ''' Sort data along axis. '''
        return np.array(sorted(data,key=lambda x:x[col]))

    def split_array(self, data_new, next_level, level, count, ind):
        ''' Split array on the node for the next level.'''
        self.book[next_level]['data_'+str(level+1)+'_'+str(count)] = data_new[:ind]
        self.book[next_level]['data_'+str(level+1)+'_'+str(count+1)] = data_new[ind:]
        return 

    def init_dict(self, data):
        self.book['level1'] = {}
        self.book['level1']['data_root'] = data

    def check_size(self, data, size):
        if len(data[:,0])<=size:
            self.nodes.append(None)

    def run(self, data, height, size=1):
        ''' Loop through levels adding nodes to array. 
            Exit either after cycling through a user defined
            number of levels (height) or when the minimum number
            of points in a cell falls below a user defined value
            (size). '''
        # Initialize first level of the dictionary
        self.init_dict(data)

        while self.level<=height:
            # Initial conditions
            count = 1
            working_level = 'level'+str(self.level)
            next_level = 'level'+str(self.level+1)
            self.book[next_level] = {}

            for file in list(self.book[working_level]):
                data = self.book[working_level][file]
                #print("{}: file {} with shape {}".format(working_level,file,np.shape(data)))
                self.check_size(data, size)
                # Sort the data, append median and find index
                data_new = self.sort_data(data,1*self.col)
                self.nodes.append(np.median(data_new[:,1*self.col],axis=0))
                ind = int(len(data_new[:,1*self.col])/2)

                # if self.level==height and np.median(data_new[:,1*self.col],axis=0)<2.4:
                #     print("------")
                #     print("{}: median = {}, n_points = {}".format(file,np.median(data_new[:,1*self.col],axis=0),len(data_new[:,1*self.col])))
                #     print(data_new[:,1*self.col])

                # Create data for next level
                self.split_array(data_new, next_level, self.level, count, ind)
                # if self.level==height and np.median(data_new[:,1*self.col],axis=0)<2.4:
                #     print(self.book[next_level]['data_'+str(self.level+1)+'_'+str(count)])
                #     print(self.book[next_level]['data_'+str(self.level+1)+'_'+str(count+1)])

                if self.level==height:
                    self.n_in_cells.append(len(self.book[next_level]['data_'+str(self.level+1)+'_'+str(count)][:,0]))
                    self.n_in_cells.append(len(self.book[next_level]['data_'+str(self.level+1)+'_'+str(count+1)][:,0]))

                count+=2

            # Update params for next loop
            self.level += 1
            self.col = np.logical_not(self.col)

            # End simulaton if single point found
            for file in list(self.book[working_level]):
                if len(self.book[working_level][file])<=size:
                    print("Single child found: Exit simulation")
                    break 
            else:
                continue
            break
        return self.nodes, self.n_in_cells


