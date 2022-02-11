from matplotlib import colors
from pylab import cm

cmap = cm.get_cmap('tab20', 20)
with open('cm_hex.txt', 'w') as file:
    for i in range(cmap.N):
        rgba = cmap(i)
        # rgb2hex accepts rgb or rgba
        file.write(colors.rgb2hex(rgba))
        file.write('\n')