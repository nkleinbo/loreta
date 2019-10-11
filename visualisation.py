from PIL import Image, ImageDraw
import math


WIDTH = 800
HEIGHT = 200
LINEWIDTH = 5
FEATUREWIDTH = 10

BLAST_FEATURE_MAPPING = {
    "Agrobacterium_circular": "agrobac",
    "Agrobacterium_plasmid": "agrobac",
    "Agrobacterium_linear": "agrobac",
    "At_Chr1": "plant",
    "At_Chr2": "plant",
    "At_Chr3": "plant",
    "At_Chr4": "plant",
    "At_Chr5": "plant",
    "Mitochondria": "mitochondria",
    "Chloroplast": "chloroplast",
    "T732"
    "LR32"
    "pAC106_tdna": "tdna",
    "pAC106_vec": "vector",
    "pAC161_tdna": "tdna",
    "pAC161_vec": "vector",
    "pGABI1_tdna": "tdna",
    "pGABI1_vec": "vector",
    "pADIS1_tdna": "tdna",
    "pADIS1_vec": "vector",
    "pROK2_vec": "vector",
    "pRok2_tdna": "tdna",
    "LB": "LB",
    "RB": "RB"
}


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



