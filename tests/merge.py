from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import easygui
font1 = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSansBold.ttf", 65)

#overlap = 100
overlap = int(input('overlap = '))
name1 = easygui.fileopenbox()
name2 = easygui.fileopenbox()
img1 = Image.open(name1)
img2 = Image.open(name2)
width1, height1 = img1.size
width2, height2 = img2.size
if height1 == height2:
    toImage = Image.new('RGBA',(width1 + width2 - overlap, height1))
    draw = ImageDraw.Draw(toImage)
    toImage.paste(img2, (width1 - overlap, 0))
    toImage.paste(img1, (0, 0))
    toImage.save('merged_hori.png')
else:
    print('heights not equal!')

def mergehori(fn_list, cuts=(0., 0.)):
    if isinstance(fn_list, list) is not True:
        print('fn_list not a list')
        return
    ok_flg = True
    imgs = []
    x_now = 0
    img = Image.open(fn_list[0])
    w, h = img.size
    h_1st = h
    target = Image.new('RGBA',(len(fn_list)*w, h))
    imgs.append(img)
    for i, fn in enumerate(fn_list[1:]):
        img = Image.open(fn)
        imgs.append(img)
        w, h = img.size
        if h != h_1st:
            print('heights not the same')
            ok_flg = False
            return
        cut_0, cut_1 = round(w*cuts[0]), round(w*cuts[1])
        target.paste(img, (x_now))        
        x_now += w - round(w*(cuts[0]+cuts[1])

    # w_total = int(np.ceil(w_total))



if __name__ = '__main__':
