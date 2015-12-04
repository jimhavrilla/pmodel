def remove_trailing_zeros(value):
    value = str(value)
    if value.find('e') != -1:
        vals = value.split('e')
        e = vals[1]
        return '{:g}'.format(float(vals[0]))+'e'+e
    vals = value.split('.')
    if (vals[0] == '0'):
        i = 0
        while vals[1][i] == '0':
            i += 1
        return '{:.{}e}'.format(float(value), len(vals[1][i:]) - 1)
    else:
        j = len(vals[0]) - 1
        while vals[0][j] == '0':
            j -= 1
        return '{:.{}e}'.format(float(value), len(vals[0][:j]))
