from PIL import Image, ImageDraw
import glob

frames = []
imgs = []
for i in range(1,155):
    imgs.append("eclipse_{}.png".format(i))

for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)

frames[0].save('movie.gif', format='GIF', append_images=frames[1:], save_all=True, duration = 100, Loop=0)