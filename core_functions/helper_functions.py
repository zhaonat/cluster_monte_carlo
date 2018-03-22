def sub2ind(array_shape, rows, cols):
    '''
    :param array_shape: (Nx,Ny)
    :param rows: x-coord
    :param cols: y-coord
    :return:
    '''
    return rows*array_shape[1] + cols

def ind2sub(array_shape, ind):
    ind[ind < 0] = -1
    ind[ind >= array_shape[0]*array_shape[1]] = -1
    rows = (ind.astype('int') / array_shape[1])
    cols = ind % array_shape[1]
    return (rows, cols)