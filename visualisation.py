from PIL import Image, ImageDraw, ImageFont
import math

WIDTH = 1200
HEIGHT = 200
OFFSET = 50
LINEWIDTH = 5
FEATUREWIDTH = 16
FONTSIZE = 12
ARROW_OFFSET = 5
FONT = "fonts/TruenoRg.otf"

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
    "T732": "adapter",
    "LR32": "adapter",
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
    "tdna": (200,0,0,1),
    "vector": (150,0,0,1),
    "adapter": (100,0,0,1),
    "plant": (0,200,0,1),
    "chloroplast": (0,155,0,1),
    "mitochondria": (0,0,255,1),
    "LB": (255,0,0,1),
    "RB": (255,0,0,1)
}

def draw_feature(draw, font, hit_dict, contig_length, top):
    x1 = math.floor((WIDTH-2*OFFSET) * hit_dict["qstart"] / contig_length)+OFFSET
    x2 = math.ceil((WIDTH-2*OFFSET) * hit_dict["qend"] / contig_length)+OFFSET
    y_line = HEIGHT / 2

    label = hit_dict["hitname"]+":\n"+str(hit_dict["sstart"])+"-"+str(hit_dict["send"])
    (textsize_x, textsize_y) = draw.multiline_textsize(label, font=font)
    x_text = int((x2+x1-textsize_x)/2)
    if top:
        y_text = y_line-FEATUREWIDTH-textsize_y
    else:
        y_text = y_line+FEATUREWIDTH
    draw.multiline_text((x_text,y_text), label, font=font, fill=COLOURS[hit_dict["f_type"]], align="center")

    if(hit_dict["sstart"] < hit_dict["send"]):
        p1 = (x1, y_line-1/2*FEATUREWIDTH)
        p2 = (x1+ARROW_OFFSET, y_line)    
        p3 = (x1, y_line+1/2*FEATUREWIDTH)
        p4 = (x2-ARROW_OFFSET, y_line+1/2*FEATUREWIDTH)
        p5 = (x2, y_line)
        p6 = (x2-ARROW_OFFSET, y_line-1/2*FEATUREWIDTH)
    else:
        p1 = (x1+ARROW_OFFSET, y_line-1/2*FEATUREWIDTH)
        p2 = (x1, y_line)    
        p3 = (x1+ARROW_OFFSET, y_line+1/2*FEATUREWIDTH)
        p4 = (x2, y_line+1/2*FEATUREWIDTH)
        p5 = (x2-ARROW_OFFSET, y_line)
        p6 = (x2, y_line-1/2*FEATUREWIDTH)
    polygon = [p1, p2, p3, p4, p5, p6]
    draw.polygon(polygon, fill=COLOURS[hit_dict["f_type"]])
        
def draw_insertion(image_name, contig_dict):
    im = Image.new("RGB", (WIDTH, HEIGHT), (255,255,255))
    contig_length = int(contig_dict["length"])
    blast_hits = contig_dict["hits"]
    draw = ImageDraw.Draw(im)
    font = ImageFont.truetype(FONT, FONTSIZE)
    y = HEIGHT / 2
    draw.line((0+OFFSET, y, WIDTH-OFFSET, y), width=LINEWIDTH, fill=(0,0,0))
    #create feature dicr:
    top = False;
    for hit in sorted(blast_hits, key = lambda i: i["qstart"]):
        hit_dict = {"hitname": hit["saccver"], 
                    "qstart": int(hit["qstart"]), 
                    "qend": int(hit["qend"]), 
                    "sstart": int(hit["sstart"]), 
                    "send": int(hit["send"]), 
                    "f_type": BLAST_FEATURE_MAPPING[hit["saccver"]]
                    }
        draw_feature(draw, font, hit_dict, contig_length, top)
        top = not top;
    im.save(image_name, "PNG")



