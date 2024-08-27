import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import json
import GeometeryClass




json_string = open('input.json').read()
json_vals = json.loads(json_string)
radius = json_vals['geometry']['cylinder_radius']
x_low_val = json_vals['plot']['x_lower_limit']
x_up_val = json_vals['plot']['x_upper_limit']


geometry = GeometeryClass.Geometery(radius)
camber, uppersurface, lower_surface = geometry.circle()


plt.figure()
plt.plot(camber[0], camber[1], label='Camber Line')
plt.plot(uppersurface[0], uppersurface[1], color='blue', label='Upper Surface')
plt.plot(lower_surface[0], lower_surface[1], color='red', label='Lower Surface')
plt.xlim(x_low_val, x_up_val)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(loc='upper right')
plt.show()

