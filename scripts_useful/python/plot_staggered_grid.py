
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Set up the 3D plot with transparent background
fig = plt.figure(figsize=(8, 6), facecolor='none')
ax = fig.add_subplot(111, projection='3d')

# Turn off coordinate grid lines
ax.grid(False)

# Remove axes
ax.set_axis_off()
# ax.set_axis_on()

# Invert the z-axis
ax.invert_zaxis()

# ax_origin=[-1.3,0.2,0]
# arrow_length = 0.8
# ax.quiver(ax_origin[0],ax_origin[1],ax_origin[2], arrow_length, 0, 0, color='k', arrow_length_ratio=0.1)  # x-axis arrow
# ax.quiver(ax_origin[0],ax_origin[1],ax_origin[2], 0, arrow_length, 0, color='k', arrow_length_ratio=0.1)  # y-axis arrow
# ax.quiver(ax_origin[0],ax_origin[1],ax_origin[2], 0, 0, arrow_length, color='k', arrow_length_ratio=0.1)  # z-axis arrow (downward)


# Node positions
nodes = {
    'V_P': [(0, 0, 0), (2, 0, 0), (0, 2, 0), (2, 2, 0), (0, 0, 2), (2, 0, 2), (0, 2, 2), (2, 2, 2)],
    'dV_dx': [(1, 0, 0), (1, 2, 0), (1, 0, 2), (1, 2, 2)],
    'dV_dy': [(0, 1, 0), (2, 1, 0), (0, 1, 2), (2, 1, 2)],
    'dV_dz': [(0, 0, 1), (2, 0, 1), (0, 2, 1), (2, 2, 1)]
}

ALPHA = 0.99    # Reduced for nodes to let halos show through
ALPHA2 = 0.6   # Increased for halos to be more visible
FSZ = 16
ScSZ1 = 400
ScSZ2 = 200
halo_size_factor = 1.8  # Increased for more visible halos

# Plot white halos on top of each node (slightly larger)
for pos in nodes['V_P']:
    ax.scatter(*pos, c='white', s=ScSZ1 * halo_size_factor, alpha=ALPHA2, edgecolor='black')
for pos in nodes['dV_dx']:
    ax.scatter(*pos, c='white', s=ScSZ2 * halo_size_factor, alpha=ALPHA2, edgecolor='black')
for pos in nodes['dV_dy']:
    ax.scatter(*pos, c='white', s=ScSZ2 * halo_size_factor, alpha=ALPHA2, edgecolor='black')
for pos in nodes['dV_dz']:
    ax.scatter(*pos, c='white', s=ScSZ2 * halo_size_factor, alpha=ALPHA2, edgecolor='black')

# Plot nodes with transparency and add labels for legend
ax.scatter(*zip(*nodes['V_P']), c='blue', s=ScSZ1, edgecolor='black', alpha=ALPHA, label=r'$V(\mathbf{x}), P(\mathbf{x},t)$')
ax.scatter(*zip(*nodes['dV_dx']), c='gray', s=ScSZ2, edgecolor='black', alpha=ALPHA, label=r'$\frac{\partial v_x(\mathbf{x},t)}{\partial x}$')
ax.scatter(*zip(*nodes['dV_dy']), c='pink', s=ScSZ2, edgecolor='black', alpha=ALPHA, label=r'$\frac{\partial v_y(\mathbf{x},t)}{\partial y}$')
ax.scatter(*zip(*nodes['dV_dz']), c='white', s=ScSZ2, edgecolor='black', alpha=ALPHA, label=r'$\frac{\partial v_z(\mathbf{x},t)}{\partial z}$')

# Add connecting lines with transparency
edges = [
    [(0,0,0), (2,0,0)], [(0,2,0), (2,2,0)], [(0,0,2), (2,0,2)], [(0,2,2), (2,2,2)],
    [(0,0,0), (0,2,0)], [(2,0,0), (2,2,0)], [(0,0,2), (0,2,2)], [(2,0,2), (2,2,2)],
    [(0,0,0), (0,0,2)], [(2,0,0), (2,0,2)], [(0,2,0), (0,2,2)], [(2,2,0), (2,2,2)]
]
for edge in edges:
    x, y, z = zip(*edge)
    ax.plot(x, y, z, 'k-', alpha=0.99)

# Add labels to key nodes with transparency
ax.text(-0.3, 0, -0.1, r'$(i,j,k)$', fontsize=FSZ, verticalalignment='bottom', alpha=ALPHA)
ax.text(2.3, -0.2, 0.01, r'$(i+1,j,k)$', fontsize=FSZ, verticalalignment='bottom', alpha=ALPHA)
ax.text(0.14, 2, -0.4, r'$(i,j+1,k)$', fontsize=FSZ, verticalalignment='top', alpha=ALPHA)
ax.text(0.2, 0, 2, r'$(i,j,k+1)$', fontsize=FSZ, verticalalignment='bottom', alpha=ALPHA)
ax.text(-0.8, 1, -0.3, r'$(i,j+1/2,k)$', fontsize=FSZ, verticalalignment='top', alpha=ALPHA)

# Add legend
legend = ax.legend(loc=(0.95, 0.35), frameon=False, fontsize=FSZ)

for i, text in enumerate(legend.get_texts()):
    if i >= 1:  # Last three labels (indices 1, 2, 3)
        text.set_fontsize(FSZ + 5)  # Increase fontsize by 4 (e.g., 16 -> 20)


# Set elevation, azimuth, and roll to match the diagram
ax.view_init(elev=17, azim=-55, roll=0)

# Save with transparent background
plt.tight_layout()
# plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.savefig('../exawave_report/staggered_grid.png', dpi=400, bbox_inches='tight', transparent=True)
plt.show()
# plt.close()

dd=1