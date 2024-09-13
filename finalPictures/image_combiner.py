import matplotlib.pyplot as plt
from PIL import Image

# Load images
img1 = Image.open('CGenFFdynamics.png')
img2 = Image.open('CGenFFPocketIndex.png')
img3 = Image.open('MorseMDCMdynamics.png')
img4 = Image.open('MorseMDCMpocketIndex.png')

# Create a figure (16x9 inches)
fig = plt.figure(figsize=(16, 9), dpi=400)

# Add first image (top-left)
ax1 = fig.add_axes([0.0, 0.48, 0.5, 0.5])  # [left, bottom, width, height]
ax1.imshow(img2)
ax1.axis('off')  # Hide axes

# Add second image (top-right)
ax2 = fig.add_axes([0.48, 0.48, 0.5, 0.5])  # [left, bottom, width, height]
ax2.imshow(img4)
ax2.axis('off')

# Add third image (bottom-left)
ax3 = fig.add_axes([0.0, 0.0, 0.5, 0.5])  # [left, bottom, width, height]
ax3.imshow(img1)
ax3.axis('off')

# Add fourth image (bottom-right)
ax4 = fig.add_axes([0.48, 0.0, 0.5, 0.5])  # [left, bottom, width, height]
ax4.imshow(img3)
ax4.axis('off')

# Save or show the result
plt.savefig('combined_image_manual.png', dpi=300)
plt.show()