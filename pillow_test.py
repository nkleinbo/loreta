from PIL import Image, ImageDraw
import sys
import math


WIDTH = 800
HEIGHT = 200
LINEWIDTH = 5
FEATUREWIDTH = 10

COLOURS = {
    "tdna": (200,0,0),
    "plant": (0,255,0),
    "LB": (255,0,0),
    "RB": (255,0,0)
}

def draw_feature(draw, feature_tupel, contig_length, f_type):
    x1 = math.floor(WIDTH * feature_tupel[0] / contig_length)
    x2 = math.ceil(WIDTH * feature_tupel[1] / contig_length)
    y = HEIGHT / 2
    draw.line((x1, y, x2, y), width=FEATUREWIDTH, fill=COLOURS[f_type])
        
def draw_insertion(image_name, contig_length, feature_dict):
    im = Image.new("RGB", (WIDTH, HEIGHT), (255,255,255))
    draw = ImageDraw.Draw(im)
    y = HEIGHT / 2
    draw.line((0, y, WIDTH, y), width=LINEWIDTH, fill=(0,0,0))
    for feature_type in feature_dict:
        feature_array = feature_dict[feature_type]
        for feature_tupel in feature_array:
            draw_feature(draw, feature_tupel, contig_length, feature_type)
    im.save(image_name, "PNG")


imagename = "test.png"
contig_l = 100000
feature_dict = {
    "tdna": [(50000,55000)],
    "LB":  [(50000, 50033)],
    "RB": [(54967, 55000)]
}

draw_insertion(imagename, contig_l, feature_dict)



