import numpy as np, math

from bokeh.io import show, output_notebook
from bokeh.models import BoxZoomTool, WheelZoomTool, Circle, HoverTool, MultiLine, Range1d, TapTool,  \
    GraphRenderer, StaticLayoutProvider, LabelSet, ColumnDataSource, \
    NodesAndAdjacentNodes, NodesAndLinkedEdges
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import figure, show
from bokeh.palettes import Spectral11

output_notebook()

def bezier(pocetak, kraj, kontrola):
    return [
        (1-s)*pocetak + 2*(1-s)*s*kontrola + s**2*kraj 
        for s in [i/100. for i in range(100)]
    ]

def nacrtaj_graf(cvorovi, ivice, putanja):
    
    putanja = list(zip(list(str(putanja[0])), putanja[1])) if putanja is not None else []
    broj_cvorova = len(cvorovi)
    oznake_cvorova = cvorovi

    plot = figure(
        title="Графовска репрезентација масеног спектра", 
        x_range=Range1d(0,150), y_range=Range1d(-2,2),
        width=1600, height=300,
        toolbar_location="right"
    )
    
    plot.axis.visible = False
    plot.xgrid.visible = False
    plot.ygrid.visible = False
        
    plot.add_tools(
        TapTool(),
        HoverTool(tooltips=None)
    )

    plot.toolbar.active_scroll = [x for x in plot.tools if isinstance(x, WheelZoomTool)][0]

    paleta_boja = Spectral11
    
    pocetni_cvorovi_ivica, krajnji_cvorovi_ivica, oznake_ivica = map(list,zip(*[(cvor, ivica[0], ivica[1]) for cvor, ivice_cvora in ivice.items() for ivica in ivice_cvora]))
    
    x = list(np.linspace(start=10, stop=140, num=broj_cvorova))
    y = [0] * broj_cvorova
    raspored_grafa = dict(zip(oznake_cvorova, zip(x, y))) 
    
    luk_ivice_sa_pozitivne_strane = True
    boje_ivica = []
    indeks_putanje = 0
    x_koordinate, y_koordinate = [], []
    visina_luka = 0
    
    for i,cvor in enumerate(cvorovi):
    
        visina_luka = 1
        iskoriscene_kombinacije_strane_luka = 0
        x1,y1 = raspored_grafa[cvor]
        
        for susedni_cvor in ivice[cvor]:
            
            # Уколико су два чвора један до другог у графу, а повезана су амино киселином
            # Између њих се исцртава ивица као права линија
            if (len(cvorovi) > (i+1)) and (susedni_cvor[0] == cvorovi[i+1]):
                visina_luka = 0
            elif visina_luka == 0:
                visina_luka = 1 * (1 if luk_ivice_sa_pozitivne_strane else -1)
            
            distanca = cvorovi.index(susedni_cvor[0]) - cvorovi.index(cvor)
            visina_luka = distanca//2 * (-1 if distanca%2 == 1 else 1)
            
            x2, y2 = raspored_grafa[susedni_cvor[0]]
            
            if indeks_putanje >= len(putanja):
                boje_ivica.append("black")
            elif cvor == putanja[indeks_putanje][1] and susedni_cvor[1] == putanja[indeks_putanje][0]:
                boje_ivica.append("green")
                indeks_putanje += 1
            else:
                boje_ivica.append("black")
            
            iskoriscene_kombinacije_strane_luka += 1
            x_koordinate.append(bezier(x1, x2, x2/2))
            y_koordinate.append(bezier(y1, y2, visina_luka))

            if visina_luka != 0:
                if iskoriscene_kombinacije_strane_luka == 2:
                    # visina_luka = (visina_luka + 1) % 4
                    iskoriscene_kombinacije_strane_luka = 0

                if luk_ivice_sa_pozitivne_strane:
                    visina_luka = -1*(visina_luka)
                    luk_ivice_sa_pozitivne_strane = False
                    iskoriscene_kombinacije_strane_luka += 1
                else:
                    visina_luka = abs(visina_luka)
                    luk_ivice_sa_pozitivne_strane = True
                    iskoriscene_kombinacije_strane_luka += 1
    
    
    graf = GraphRenderer()
    
    graf.layout_provider = StaticLayoutProvider(graph_layout=raspored_grafa)

    graf.node_renderer.data_source.data = {
        "index": oznake_cvorova,
        "boja": list(paleta_boja * math.ceil(broj_cvorova/len(paleta_boja)))[:broj_cvorova]
    }
    
    graf.node_renderer.glyph = Circle(radius=1.2, fill_color="boja")
    graf.node_renderer.selection_glyph = Circle(radius=1.2, fill_color="white")
    graf.node_renderer.hover_glyph = Circle(radius=1.2, fill_color="white")

    graf.edge_renderer.data_source.data = {
        'start':pocetni_cvorovi_ivica,
        'end':krajnji_cvorovi_ivica,
        'xs': x_koordinate,
        'ys': y_koordinate,
        'boje_ivica': boje_ivica,
        'sirina_linije':[2 if x == "green" else 1 for x in boje_ivica]
    }
    
    graf.edge_renderer.glyph = MultiLine(line_color="boje_ivica", line_alpha=0.8, line_width='sirina_linije')
    graf.edge_renderer.selection_glyph = MultiLine(line_color="boje_ivica", line_width=4)

    graf.selection_policy = NodesAndLinkedEdges()
    graf.inspection_policy = NodesAndAdjacentNodes()
    
    x,y = zip(*raspored_grafa.values())
    y = [-0.04 for i in y]
    
    # Додаје ознаке чворова
    cvorovi_source = ColumnDataSource({'x':list(x),'y':y, 'text':list(map(str, cvorovi))})
    labele_cvorova = LabelSet(x="x", y="y", text="text", source=cvorovi_source, text_font_size="7pt", text_color="black", text_align="center")
    
    # Додаје ознаке ивица 
    ivice_source = ColumnDataSource(data=dict(x=[x[49] for x in x_koordinate], y=[z if z > 0 else z - 0.2 for z in [y[49] for y in y_koordinate]], text=oznake_ivica))
    labele_ivica = LabelSet(x="x", y="y", text="text", source=ivice_source, text_font_size="8pt", text_color="black", text_align="center")
    
    plot.renderers.append(graf)
    plot.renderers.append(labele_cvorova)
    plot.renderers.append(labele_ivica)

    x_callback = CustomJS(args=dict(labels=labele_cvorova, init_font_size=labele_cvorova.text_font_size[:-2], init_xrange=140), code="""
        let xzoom = (init_font_size * init_xrange) / (cb_obj.end - cb_obj.start);
        labels.text_font_size = String(xzoom) + 'pt';
    """)
    plot.x_range.js_on_change('start', x_callback)


    x_callback = CustomJS(args=dict(labels=labele_ivica, init_font_size=labele_ivica.text_font_size[:-2], init_xrange=140), code="""
        let xzoom = (init_font_size * init_xrange) / (cb_obj.end - cb_obj.start);
        labels.text_font_size = String(xzoom) + 'pt';
    """)
    plot.x_range.js_on_change('start', x_callback)
    
    # Уклања BoxZoomTool
    for tool in [tool for tool in plot.tools if isinstance(tool, BoxZoomTool)]:
        plot.remove_tools(tool)

    # Уклања уколико постоји било који дуплирани алат
    tool_classes = [x.__class__ for x in plot.tools]
    for x in set(tool_classes):
        if tool_classes.count(x) > 1:
            plot.remove_tools(plot.tools[tool_classes.index(x)])
        
    show(plot)