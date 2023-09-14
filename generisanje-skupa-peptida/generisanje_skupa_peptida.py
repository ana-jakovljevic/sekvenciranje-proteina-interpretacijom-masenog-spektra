from utils import *
import random
import gc

verovatnoce_jona_po_rangovima = pd.read_csv("verovatnoce-jona-po-rangovima.csv", index_col="ion").drop("y-h2o")
celobrojni_pomeraj_fragmentnih_jona = {
    "y": 19,
    "b": 1, 
    "b-h2o": 17,
    #"y-h2o": 1,
    "a": -27,
    "b-nh3": -16,
    "y-nh3": 2,
    "b-h2o-h2o": -35,
}


def generisi_slucajan_peptid(duzina=None):
    if duzina is None:
        duzina = random.randint(8,12)
        
    return "".join([random.choice(list(amino_kiseline.keys())) for i in range(duzina)])

def generisi_fragmentne_jone_peptida_deprecated(peptid, izabrani_joni, skaliraj_rangove):
    fragmentni_joni_peptida = []
    
    for i in range(len(peptid)-1):
        prefiksna_masa_peptida = Peptid(peptid[0:i]).masa
        sufiksna_masa_peptida = Peptid(peptid[i:]).masa  
    
        # sa kojom se verovatnocom jon javlja => nece se javiti svaki put, vec samo ponekad
        # a mi pokusavamo da simuliramo to ponekad
        # za svaki prefiks/sufiks => za svaki tip jona generisi ga sa odgovarajucom verovatnocom => ali je potrebno izabrati u kom opsegu intenziteta se nalazi ta verovatnoca koju zelimo da generisemo (dakle imamo dva parametra)
        for rang, verovatnoce_jona in enumerate(verovatnoce_jona_po_rangovima.to_dict().items()):
            izabrani_jon = random.choices(list(verovatnoce_jona[1].keys()), verovatnoce_jona[1].values(), k=1)[0]            
            
            if skaliraj_rangove:
                intenzitet = random.randint(10 + (40 - rang*10) * 250, 10 + (40 - rang*10 + 1 ) * 250)
            else:
                intenzitet = random.randint(10 + (3 - rang) * 250, 10 + (3 - rang + 1 ) * 250)
                
            if izabrani_jon == "noise":
                fragmentni_joni_peptida.append(generisi_sum(intenzitet = intenzitet))
            else:
                pomeraj = celobrojni_pomeraj_fragmentnih_jona[izabrani_jon]

                # sufiksni peptidi
                if np.any([izabrani_jon.startswith(oznaka) for oznaka in ["x","y","z"]]):
                    if (sufiksna_masa_peptida + pomeraj) >= 91 and (sufiksna_masa_peptida + pomeraj) <= 2000:
                        fragmentni_joni_peptida.append([sufiksna_masa_peptida + pomeraj, intenzitet])

                # prefiksni peptidi
                elif np.any(izabrani_jon.startswith(oznaka) for oznaka in ["a","b","c"]):
                    if (prefiksna_masa_peptida + pomeraj) >= 91 and (prefiksna_masa_peptida + pomeraj) <= 2000:
                        fragmentni_joni_peptida.append([prefiksna_masa_peptida + pomeraj, intenzitet])
                        
    return fragmentni_joni_peptida

def generisi_fragmentne_jone_peptida(peptid, izabrani_joni, skaliraj_rangove):
    izabrane_m_z = []
    fragmentni_joni_peptida = []
    pozicije_prekida = list(range(len(peptid)))

    for rang, verovatnoce_jona in verovatnoce_jona_po_rangovima.to_dict().items():

        i = 0
        while i < 10:
            
            intenzitet = random.randint((4 - int(rang)) * 500, (4 - int(rang) + 1 ) * 500)

            izabrani_jon = random.choices(list(verovatnoce_jona.keys()), verovatnoce_jona.values(), k=1)[0]
            izabrana_pozicija = random.choice(pozicije_prekida)

            pomeraj = celobrojni_pomeraj_fragmentnih_jona[izabrani_jon] if izabrani_jon != "noise" \
                else generisi_sum()

            # prefiks
            if np.any([izabrani_jon.startswith(oznaka) for oznaka in ["x","y","z"]]):
                masa_fragmenta = Peptid(peptid[izabrana_pozicija:]).masa 
            # sufiks
            elif np.any(izabrani_jon.startswith(oznaka) for oznaka in ["a","b","c"]): 
                masa_fragmenta = Peptid(peptid[izabrana_pozicija:]).masa 
                masa_fragmenta = Peptid(peptid[:izabrana_pozicija]).masa
            # sum
            else: 
                masa_fragmenta = 0

            m_z = masa_fragmenta + pomeraj

            if m_z in izabrane_m_z:
                continue
            else:
                i += 1 

            if (masa_fragmenta + pomeraj) >= 91 and (masa_fragmenta + pomeraj) <= 2000:
                fragmentni_joni_peptida.append([m_z, intenzitet])
                
            izabrane_m_z.append(m_z)
            
    return fragmentni_joni_peptida
            
def generisi_sum():

    # simulacija kontaminanta (sum)    
    peptid = generisi_slucajan_peptid()
    
    for i in range(len(peptid)):
        prefiksna_masa_peptida = Peptid(peptid[0:i]).masa
        sufiksna_masa_peptida = Peptid(peptid[i:]).masa  
            
        m_z = random.choice([prefiksna_masa_peptida, sufiksna_masa_peptida])
        
    return m_z

def generisi_trening_skup(broj_peptida, skaliraj_rangove=False):
    izabrani_joni = []
    peptidi = [generisi_slucajan_peptid() for i in range(broj_peptida)]
   
    trening_skup = [
        (Peptid(peptid), MaseniSpektar(
            pd.DataFrame(generisi_fragmentne_jone_peptida(peptid, izabrani_joni, skaliraj_rangove), columns=["m_z","intenzitet"])))
     for peptid in peptidi
    ]
    
    return izabrani_joni, TreningSkup(trening_skup)

def dodeli_rangove_spektru_dancik(spektar, roditeljska_masa):
    
    spektar = spektar.sort_values(by="intenzitet", ascending=False).reset_index(drop=True)

    K = roditeljska_masa // 100 # velicina bina 
    broj_rangova = len(spektar) // K + 1 

    for j in range(0, broj_rangova):
        spektar.loc[j*K:(j+1)*K, "rang"] = j + 1
        
    spektar.rang = spektar.rang.astype('int')
    return spektar

def generisi_eksperimentalnu_verovatnoca_jona_po_rangovima_dancik(trening_skup):
    objasnjeni_pikovi_ranga = {}
    neobjasnjeni_pikovi_ranga = {}
    broj_odgovarajucih = {jon:{} for jon in celobrojni_pomeraj_fragmentnih_jona.keys()} 

    for i in range(0,len(trening_skup.peptidi_i_spektri)):    
        peptid = trening_skup.peptidi_i_spektri[i][0] #.peptid
        spektar = trening_skup.peptidi_i_spektri[i][1].mase_intenziteti

        spektar = dodeli_rangove_spektru_dancik(spektar, peptid.masa)
        peptid = peptid.peptid
        duzina_peptida = len(peptid)

        for j in range(duzina_peptida):
            prefiksna_masa = Peptid(peptid[0:j]).masa
            sufiksna_masa = Peptid(peptid[j:]).masa 

            for jon, pomeraj in celobrojni_pomeraj_fragmentnih_jona.items():
                for rang in range(1, max(spektar.rang)+1):

                    if rang not in broj_odgovarajucih[jon]:
                        broj_odgovarajucih[jon][rang] = 0

                    # sufiksni peptidi
                    if np.any([jon.startswith(oznaka) for oznaka in ["x","y","z"]]):
                        # broj_pozicija[jon][rang] += (1 if ((sufiksna_masa + pomeraj >= 91) and (sufiksna_masa + pomeraj <= 2000)) else 0)
                        sufiksi_u_spektru = spektar[
                            (spektar["m_z"] >= sufiksna_masa + pomeraj - 0.5) & 
                            (spektar["m_z"] <= sufiksna_masa + pomeraj + 0.5) &
                            (spektar["rang"] == rang)
                        ]

                        broj_odgovarajucih[jon][rang] += len(sufiksi_u_spektru)
                        objasnjeni_pikovi_ranga[rang] = objasnjeni_pikovi_ranga.get(rang, []) + list(sufiksi_u_spektru.m_z.values)

                    # prefiksni peptidi
                    elif np.any(jon.startswith(oznaka) for oznaka in ["a","b","c"]):
                        # broj_pozicija[jon][rang] += (1 if ((prefiksna_masa + pomeraj >= 91) and (prefiksna_masa + pomeraj <= 2000)) else 0)
                        prefiksi_u_spektru = spektar[
                            (spektar["m_z"] >= prefiksna_masa + pomeraj - 0.5) & 
                            (spektar["m_z"] <= prefiksna_masa + pomeraj + 0.5) &
                            (spektar["rang"] == rang)
                        ]

                        broj_odgovarajucih[jon][rang] += len(prefiksi_u_spektru)
                        objasnjeni_pikovi_ranga[rang] = objasnjeni_pikovi_ranga.get(rang, []) + list(prefiksi_u_spektru.m_z.values)

        for rang in range(1, max(spektar.rang)+1):
            neobjasnjeni_pikovi_ranga[rang] = neobjasnjeni_pikovi_ranga.get(rang, 0) + len(list(set(spektar.loc[spektar["rang"] == rang, "m_z"].values).difference(set(objasnjeni_pikovi_ranga[rang]))))
            
    eksperimentalna_verovatnoca_jona_po_rangovima = {
        jon: {
            rang: round(broj_odgovarajucih[jon][rang]/len(objasnjeni_pikovi_ranga[rang]), 5) for rang in po_rangovima
        } for jon, po_rangovima in broj_odgovarajucih.items()
    }
    
    broj_rangova = max(pd.DataFrame(eksperimentalna_verovatnoca_jona_po_rangovima).index)
    
    eksperimentalna_verovatnoca_jona_po_rangovima["noise"] = {
        rang: round(neobjasnjeni_pikovi_ranga[rang]/len(objasnjeni_pikovi_ranga[rang]), 5) for rang in range(1, broj_rangova+1)
    }
    
    return eksperimentalna_verovatnoca_jona_po_rangovima

def dodeli_rangove_spektru_generating(spektar):
    
    spektar = spektar.sort_values(by="intenzitet", ascending=False).reset_index(drop=True)
    spektar["rang"] = spektar.index + 1 
    spektar.rang = spektar.rang.astype('int')
    spektar.set_index(keys=["rang"], drop=True, inplace=True)
    
    return spektar

def generisi_eksperimentalnu_verovatnoca_jona_po_rangovima_generating(trening_skup):
    broj_pikova_po_rangu = {}
    broj_odgovarajucih = {}
    objasnjeni_pikovi_ranga = {}
    neobjasnjeni_pikovi_ranga = {}

    for i in range(0, len(trening_skup.peptidi_i_spektri)):   
        peptid = trening_skup.peptidi_i_spektri[i][0].peptid
        spektar = dodeli_rangove_spektru_generating(trening_skup.peptidi_i_spektri[i][1].mase_intenziteti)
        
        duzina_peptida = len(peptid)
        broj_rangova = max(spektar.index.get_level_values(level="rang")) + 1
        
        for rang in range(1,  broj_rangova+1):
            if rang not in spektar.index.get_level_values(level="rang"):
                continue

            #broj_pikova_po_rangu[rang] += len(spektar.loc[[rang]])
            broj_pikova_po_rangu[rang] = broj_pikova_po_rangu.get(rang, 0) + len(spektar.loc[[rang]])

        for j in range(duzina_peptida):
            prefiksna_masa = Peptid(peptid[0:j]).masa
            sufiksna_masa = Peptid(peptid[j:]).masa 

            for rang in range(1, broj_rangova):    
                for jon, pomeraj in celobrojni_pomeraj_fragmentnih_jona.items():

                    if rang not in spektar.index.get_level_values(level="rang"):
                        continue

                    if rang not in broj_odgovarajucih:
                        broj_odgovarajucih[rang] = {}

                    if jon not in broj_odgovarajucih[rang]:
                        broj_odgovarajucih[rang][jon] = 0

                    # sufiksni peptidi
                    if np.any([jon.startswith(oznaka) for oznaka in ["x","y","z"]]):
                        # broj_pozicija[jon][rang] += (1 if ((sufiksna_masa + pomeraj >= 91) and (sufiksna_masa + pomeraj <= 2000)) else 0)
                        sufiksi_u_spektru = spektar.loc[[rang]][lambda df: (df.m_z >= sufiksna_masa + pomeraj - 0.5) & (df.m_z <= sufiksna_masa + pomeraj + 0.5)]
                        broj_odgovarajucih[rang][jon] += len(sufiksi_u_spektru)
                        objasnjeni_pikovi_ranga[rang] = objasnjeni_pikovi_ranga.get(rang, []) + list(sufiksi_u_spektru.m_z.values)

                    # prefiksni peptidi
                    elif np.any(jon.startswith(oznaka) for oznaka in ["a","b","c"]):
                        # broj_pozicija[jon][rang] += (1 if ((prefiksna_masa + pomeraj >= 91) and (prefiksna_masa + pomeraj <= 2000)) else 0)
                        prefiksi_u_spektru = spektar.loc[[rang]][lambda df: (df.m_z >= prefiksna_masa + pomeraj - 0.5) & (df.m_z <= prefiksna_masa + pomeraj + 0.5)]
                        broj_odgovarajucih[rang][jon] += len(prefiksi_u_spektru)  
                        objasnjeni_pikovi_ranga[rang] = objasnjeni_pikovi_ranga.get(rang, []) + list(prefiksi_u_spektru.m_z.values)

        for rang in range(1, max(spektar.index)+1):
            neobjasnjeni_pikovi_ranga[rang] = neobjasnjeni_pikovi_ranga.get(rang, 0) + len(list(set(spektar.loc[[rang]]["m_z"].values).difference(set(objasnjeni_pikovi_ranga[rang]))))

    eksperimentalna_verovatnoca_jona_po_rangovima = {
        rang: {
            jon: broj_odgovarajucih[rang][jon]/len(objasnjeni_pikovi_ranga[rang]) for jon in po_jonima
        } for rang, po_jonima in broj_odgovarajucih.items()
    }

    broj_rangova = max(pd.DataFrame(eksperimentalna_verovatnoca_jona_po_rangovima).index)
    
    eksperimentalna_verovatnoca_jona_po_rangovima["noise"] = {
        rang: round(neobjasnjeni_pikovi_ranga[rang]/len(objasnjeni_pikovi_ranga[rang]), 5) for rang in range(1, broj_rangova+1)
    }

    #for i in range(1, 1 + len(eksperimentalna_verovatnoca_jona_po_rangovima)):
    #    eksperimentalna_verovatnoca_jona_po_rangovima[i]["noise"] = 1.0 - sum(eksperimentalna_verovatnoca_jona_po_rangovima[i].values())
    #    eksperimentalna_verovatnoca_jona_po_rangovima[i]["noise"] = max(0, eksperimentalna_verovatnoca_jona_po_rangovima[i]["noise"])
    
    skor_po_rangovima = {
    rang: {
            jon: math.log(eksperimentalna_verovatnoca_jona_po_rangovima[rang][jon]/eksperimentalna_verovatnoca_jona_po_rangovima[rang]["noise"]) \
            if (eksperimentalna_verovatnoca_jona_po_rangovima[rang][jon] != 0) and (eksperimentalna_verovatnoca_jona_po_rangovima[rang]["noise"] != 0) \
            else math.log(eksperimentalna_verovatnoca_jona_po_rangovima[rang][jon]) if (eksperimentalna_verovatnoca_jona_po_rangovima[rang][jon] != 0) \
            else 0
                for jon in po_jonima
        } for rang, po_jonima in eksperimentalna_verovatnoca_jona_po_rangovima.items()
    }

    skor_po_rangovima_po_jonima = pd.DataFrame.from_dict(skor_po_rangovima, orient='index')

    return eksperimentalna_verovatnoca_jona_po_rangovima, skor_po_rangovima_po_jonima

class Peptid:
    def __init__(self, peptid):
        self.peptid = peptid
        self.masa = int(sum([amino_kiseline[self.peptid[i]] for i in range(len(self.peptid))]))
        
        self.prefiksi = [self.peptid[:i+1] for i in range(len(self.peptid))]
        self.sufiksi = [self.peptid[i:] for i in range(len(self.peptid))]
        
        self.vektor_prefiksnih_masa = list(np.cumsum([0] + [amino_kiseline[self.peptid[i]] for i in range(len(self.peptid))]))
        self.vektor_sufiksnih_masa = list(np.cumsum([0] + [amino_kiseline[self.peptid[-(i+1)]] for i in range(len(self.peptid))]))
        
        self.peptidni_vektor = [1 if i in self.vektor_prefiksnih_masa else 0 for i in range((int( self.vektor_prefiksnih_masa[-1]) + 1))]
                        
    def __str__(self):
        return self.peptid
    
    def __repr__(self):
        return self.peptid
    
    def __add__(self, drugi_peptid):
        return Peptid(str(self.peptid) + str(drugi_peptid))
    
    def __lt__(self, drugi_peptid):
        return self.peptid < drugi_peptid.peptid
    
    def __le__(self, drugi_peptid):
        return self.peptid <= drugi_peptid.peptid
    
    def __eq__(self, drugi_peptid):
        return self.peptid == drugi_peptid.peptid
    
    def __ne__(self, drugi_peptid):
        return self.peptid != drugi_peptid.peptid
    
    def __gt__(self, drugi_peptid):
        return self.peptid > drugi_peptid.peptid
    
    def __ge__(self, drugi_peptid):
        return self.peptid >= drugi_peptid.peptid
    
    def __hash__(self):
        return self.peptid.__hash__()
    
    def __len__(self):
        return len(self.peptid)
    
    def peptidni_vektor_u_peptid(clf, vektor):
        prefiksne_mase = [i+1 for i, masa in enumerate(vektor) if masa]
        sekvenca_amino_kiselina = prefiksne_mase - np.concatenate(([0], prefiksne_mase[:-1]))
        
class MaseniSpektar:
    # Uvodimo intenzitete u pricu pa spektar ne sadrzi samo mase, i sada se predstavlja dvodimenzionom strukturom (izabran je DataFrame) koji sadrzi dve kolone: m_z i intentizteti
    
    def __init__(self, spektar):
        self.mase_intenziteti = spektar.sort_values(by="intenzitet", ascending=False).reset_index(drop=True)
        self.m_z = sorted(list(spektar["m_z"].astype('int64')))
        self.intenziteti = list(spektar["intenzitet"])        
    
    def nacrtaj(self):
        fig, ax = plt.subplots(figsize=(15, 5))

        ax.bar(self.m_z, height=self.intenziteti, width=2, linewidth=2)
        #oznake_pikova = [f"{self.m_z[i]}: {self.intenziteti[i]}" if self.intenziteti[i] >= 250 else 0 for i in range(len(self.intenziteti))]

        #for indeks, oznaka in enumerate(list(oznake_pikova)):
        #    if oznaka != 0:
        #        ax.annotate(oznaka, (self.m_z[indeks] - 10, self.intenziteti[indeks] + 20))

        plt.show()
        
class TreningSkup:

    def __init__(self, peptidi_i_spektri):
        self.peptidi_i_spektri = peptidi_i_spektri
    
    @classmethod
    def izracunaj_pomeraje(cls, peptidne_mase, spektralne_mase, naelektrisanje):
        return [(p_i - naelektrisanje * s_j) for s_j in spektralne_mase for p_i in peptidne_mase]
    
    def funkcija_ucestalosti_pomeraja(self, naelektrisanja, prefiks):
        fig, axs = plt.subplots(len(naelektrisanja), 1, figsize=(50, len(naelektrisanja)*10))
       
        axs = axs if len(naelektrisanja) != 1 else [axs]
       
        for i, naelektrisanje in enumerate(naelektrisanja):
            pomeraji = []

            for (peptid, spektar) in self.peptidi_i_spektri:
                peptidne_mase = peptid.vektor_prefiksnih_masa if prefiks else peptid.vektor_sufiksnih_masa
                pomeraji += TreningSkup.izracunaj_pomeraje(peptidne_mase, spektar.m_z, naelektrisanje)

            axs[i].hist(pomeraji, bins=np.arange(-50, 50, 1), density=True, stacked=True, align='left')
            axs[i].set_xticks(np.arange(-50, 50, 1))
            axs[i].set_title(f"naelektrisanje={naelektrisanje}; {'prefiksi' if prefiks else 'sufiksi'}")

        plt.show()        
   
    def funkcija_ucestalosti_pomeraja_po_grupama_intenziteta(self, naelektrisanje, prefiks):   
        
        pomeraji = dict({})

        for i, (peptid, spektar) in enumerate(trening_skup.peptidi_i_spektri):
            # grupe = dict({})

            sortirani_spektar = spektar.mase_intenziteti.sort_values(by="intenzitet", ascending=False)
            peptidne_mase = peptid.vektor_prefiksnih_masa if prefiks else peptid.vektor_sufiksnih_masa

            K = peptid.masa // 100  # max(spektar.m_z) => TODO: roditeljska masa spektra ? 
            broj_rangova = len(sortirani_spektar) // K + 1 

            for rang in range(0, broj_rangova):
                grupa = sortirani_spektar.iloc[rang*K:(rang+1)*K]
                pomeraji[rang] = pomeraji.get(rang, []) + TreningSkup.izracunaj_pomeraje(peptidne_mase, list(grupa.m_z.values), naelektrisanje)
           
        broj_rangova = len(pomeraji)  # + 1 
        fig, axs = plt.subplots(broj_rangova, 1, figsize=(50, 15*broj_rangova))

        for k in range(broj_rangova):
            pomeraji_manji_od_k = [x for i in range(0, k) for x in pomeraji.get(i, [])]
            pomeraji_veci_jednaki_od_k = [x for i in range(k, broj_rangova) for x in pomeraji.get(i,[])]

            histogram_manjih = np.histogram(pomeraji_manji_od_k, bins=np.arange(-40, 40, 1))
            histogram_vecih = np.histogram(pomeraji_veci_jednaki_od_k, bins=np.arange(-40, 40, 1))

            const = 0.5
            x = np.concatenate((const*(histogram_manjih[1][1:] + histogram_manjih[1][:-1]), const*(histogram_vecih[1][1:] + histogram_vecih[1][:-1])))
            x = [i-const for i in x]
            y = np.concatenate((histogram_manjih[0] * -1, histogram_vecih[0]))

            broj_oznaka = len(histogram_manjih[1])
            axs[k].bar(x=x, height=y, width=1, align='center')
            
            oznake_pikova = [
                mase_fragmentnih_jona_prefiksa.get(-1 * int(i), '') if prefiks else mase_fragmentnih_jona_sufiksa.get(-1 * int(i), '')\
                for i in x[broj_oznaka:]
            ]

            for indeks, oznaka in enumerate(list(oznake_pikova)):
                if oznaka != '':
                    axs[k].annotate(oznaka, (x[broj_oznaka+indeks], y[broj_oznaka+indeks] + 100), fontsize='large', ha='center')

            axs[k].set_xticks(np.arange(-41, 41, 1))
            axs[k].set_yticks(np.arange(
                  int(np.ceil((-1 * max([abs(i) for i in y]) - 1000)/1000)*1000), 
                  max([abs(i) for i in y]) + 1000,
                  1000
                )
             )
            axs[k].set_ylim(-1 * max([abs(i) for i in y]) - 1000, max([abs(i) for i in y]) + 1000)
            axs[k].axhline(y=0, xmin=0, xmax=1, color="black")
            axs[k].set_title(f"rang >= {k+1}, rang < {k+1}")

        plt.show()