"""
build_vogeltrek.py — GeoPandas-pipeline voor vogeltrek.html (kraanvogel, Grus grus)

Gebruik:
    python build_vogeltrek.py <map_met_natural_earth_geojson> <uitvoermap>

Invoer (Natural Earth 1:50m, public domain):
    ne_50m_admin_0_countries.geojson, ne_50m_lakes.geojson,
    ne_50m_rivers_lake_centerlines.geojson

Uitvoer (GeoJSON / JSON, WGS84):
    landen.geojson     – landen geclipt op Europa/N-Afrika/Midden-Oosten, vereenvoudigd
    water.geojson      – meren + rivieren (subtiele achtergrond)
    rustplaatsen.geojson – broed-, rust- en overwinteringsgebieden met piektellingen
    flyways.geojson    – flyway-hartlijnen (Chaikin-afgevlakt) + lengte in km (EPSG:3035)
    corridors.geojson  – 150 km gebufferde corridors (EPSG:3035 -> 4326)
    tracks.json        – gemodelleerde individuen (pad + tijdstempels in dagen sinds 1 aug)
    bezetting.json     – gemodelleerde bezetting per rustplaats per dag + aandeel 'in trek'

LET OP: de tracks zijn een SIMULATIE op basis van bekende flyways, rustplaatsen en
fenologie (literatuur / tellingen), geen echte GPS-zenderdata. De echte waarnemingen
(GBIF / eBird) worden live in de browser opgehaald.
"""
import json
import math
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point, box

RNG = np.random.default_rng(20261031)
BBOX = box(-25, 4, 60, 72)          # lon/lat-uitsnede kaart
LAEA = "EPSG:3035"                    # ETRS89-LAEA Europa, voor afstanden/buffers
N_BIRDS = 640

# ---------------------------------------------------------------------------
# 1. Gebieden (benaderde centroiden; piek = orde van grootte max. dagtelling)
# ---------------------------------------------------------------------------
SITES = {
    # broedgebieden
    "b_scan": dict(naam="Broedgebied Zweden/Noorwegen", land="SE/NO", lon=15.5, lat=61.8, rol="broed", piek=None, sx=3.0, sy=2.2,
                   info="Moerassen en veenbossen in Scandinavië; grootste deel van de West-Europese populatie."),
    "b_fin": dict(naam="Broedgebied Finland", land="FI", lon=26.5, lat=63.2, rol="broed", piek=None, sx=2.6, sy=1.8,
                  info="Uitgestrekte venen; vogels trekken grotendeels via Estland en Hongarije."),
    "b_balt": dict(naam="Broedgebied Baltische staten/Polen", land="EE/LV/LT/PL", lon=24.0, lat=56.3, rol="broed", piek=None, sx=2.4, sy=1.6,
                   info="Dichte broedpopulatie in hoogvenen en moerasbossen."),
    "b_ger": dict(naam="Broedgebied Noordoost-Duitsland", land="DE", lon=13.2, lat=53.4, rol="broed", piek=None, sx=1.6, sy=0.8,
                  info="Mecklenburg-Vorpommern en Brandenburg; sterk toegenomen populatie."),
    "b_nl": dict(naam="Broedgebied Noord-Nederland", land="NL", lon=6.5, lat=52.9, rol="broed", piek=None, sx=0.5, sy=0.3,
                 info="Sinds 2001 weer broedend (o.a. Fochteloërveen, Dwingelderveld, Engbertsdijksvenen)."),
    "b_rus": dict(naam="Broedgebied West-Rusland/Wit-Rusland", land="RU/BY", lon=34.0, lat=57.0, rol="broed", piek=None, sx=4.5, sy=2.6,
                  info="Oostelijke populatie; trekt via Oekraïne en Turkije naar Israël en Ethiopië."),
    # rustplaatsen
    "hornborga": dict(naam="Hornborgasjön", land="SE", lon=13.55, lat=58.32, rol="rust", piek=25000, sx=0.15, sy=0.1,
                      info="Beroemde voorjaarsrustplaats ('kraanvogeldans'), piek eind maart – half april."),
    "matsalu": dict(naam="Matsalu", land="EE", lon=23.7, lat=58.75, rol="rust", piek=25000, sx=0.25, sy=0.15,
                    info="Verzamelplaats aan de Estse kust, piek september."),
    "rugen": dict(naam="Rügen-Bock-regio", land="DE", lon=13.05, lat=54.45, rol="rust", piek=70000, sx=0.25, sy=0.12,
                  info="Grootste verzamelplaats aan de Oostzee, piek eind september – oktober."),
    "linum": dict(naam="Linum / Rhinluch", land="DE", lon=12.87, lat=52.76, rol="rust", piek=80000, sx=0.15, sy=0.1,
                  info="Brandenburg; tienduizenden slaapplaatsvogels in oktober."),
    "diepholz": dict(naam="Diepholzer Moorniederung", land="DE", lon=8.6, lat=52.55, rol="rust", piek=120000, sx=0.2, sy=0.12,
                     info="Belangrijkste rustgebied van Noordwest-Duitsland, piek oktober – november; ook in het voorjaar."),
    "fochtelo": dict(naam="Fochteloërveen / Drenthe", land="NL", lon=6.4, lat=53.0, rol="rust", piek=6000, sx=0.15, sy=0.1,
                     info="Nederlands broed- en pleisterplaats; doortrek in okt–nov en feb–mrt goed zichtbaar boven Oost-Nederland."),
    "der": dict(naam="Lac du Der-Chantecoq", land="FR", lon=4.75, lat=48.58, rol="rust", piek=200000, sx=0.2, sy=0.12,
                info="Grootste stuwmeer van Frankrijk; tot ruim 200.000 kraanvogels in november, ook in feb–mrt."),
    "arjuzanx": dict(naam="Arjuzanx (Landes)", land="FR", lon=-0.85, lat=44.12, rol="winter", piek=50000, sx=0.3, sy=0.2,
                     info="Overwinteringsgebied in de Landes; steeds meer vogels blijven hier de winter over."),
    "gallocanta": dict(naam="Laguna de Gallocanta", land="ES", lon=-1.5, lat=40.97, rol="rust", piek=100000, sx=0.15, sy=0.1,
                       info="Steppemeer in Aragón; doortrekpiek in november en februari."),
    "extremadura": dict(naam="Dehesa's van Extremadura", land="ES", lon=-6.0, lat=39.2, rol="winter", piek=130000, sx=1.0, sy=0.6,
                        info="Eikenweides (eikels!) — belangrijkste overwinteringsgebied van de westelijke flyway."),
    "biebrza": dict(naam="Biebrza-moerassen", land="PL", lon=22.6, lat=53.5, rol="rust", piek=15000, sx=0.3, sy=0.15,
                    info="Grote laagveenmoerassen in Noordoost-Polen."),
    "hortobagy": dict(naam="Hortobágy", land="HU", lon=21.1, lat=47.58, rol="rust", piek=150000, sx=0.3, sy=0.15,
                      info="Poesta en visvijvers; recordaantallen (100–150 duizend) in oktober – november."),
    "ichkeul": dict(naam="Ichkeul / Noord-Tunesië", land="TN", lon=9.65, lat=37.1, rol="winter", piek=20000, sx=0.6, sy=0.3,
                    info="Overwintering in Noord-Afrika (Tunesië/Algerije)."),
    "algerije": dict(naam="Hoogvlakten Noord-Algerije", land="DZ", lon=5.5, lat=35.7, rol="winter", piek=15000, sx=1.0, sy=0.4,
                     info="Zoutvlakten (chotts) en akkers; deels ook Libië."),
    "polesia": dict(naam="Polesië (Pripjat-moerassen)", land="BY/UA", lon=27.5, lat=51.8, rol="rust", piek=30000, sx=1.0, sy=0.4,
                    info="Enorm moeras- en riviergebied op de grens van Wit-Rusland en Oekraïne."),
    "askania": dict(naam="Askania-Nova / Sivasj", land="UA", lon=33.9, lat=46.4, rol="rust", piek=60000, sx=0.4, sy=0.2,
                    info="Steppe en lagunes aan de Zwarte Zee, piek september – oktober."),
    "anatolie": dict(naam="Centraal-Anatolië (Sultan Sazlığı)", land="TR", lon=35.25, lat=38.3, rol="rust", piek=20000, sx=0.5, sy=0.3,
                     info="Steppemeren op de route naar het Midden-Oosten."),
    "hula": dict(naam="Hula-vallei", land="IL", lon=35.6, lat=33.1, rol="winter", piek=45000, sx=0.08, sy=0.06,
                 info="Agamon Hula: tot ~45.000 overwinteraars, met grootschalige bijvoedering."),
    "ethiopie": dict(naam="Ethiopische Riftvallei", land="ET", lon=38.9, lat=8.6, rol="winter", piek=20000, sx=1.0, sy=0.8,
                     info="Zuidelijkste overwinteringsgebied van de oostelijke populatie."),
}

# Flyways: 0 = West-Europees, 1 = Baltisch-Hongaars, 2 = Oost-Europees/Midden-Oosten
FLYWAYS = [
    dict(id=0, naam="West-Europese flyway", kleur="#4fd1c5",
         hartlijn=["b_scan", "rugen", "diepholz", "der", "gallocanta", "extremadura"]),
    dict(id=1, naam="Baltisch-Hongaarse flyway", kleur="#f6ad55",
         hartlijn=["b_fin", "matsalu", "biebrza", "hortobagy", "ichkeul"]),
    dict(id=2, naam="Oost-Europese flyway", kleur="#f472b6",
         hartlijn=["b_rus", "polesia", "askania", "anatolie", "hula", "ethiopie"]),
]

# Groepen: (aandeel, flyway, najaarsroute, voorjaarsroute, vertrek-dag (sinds 1 aug), terugkeer-vertrek-dag)
# Verblijfsduur (dagen) per rustplaats: (gemiddeld, sd)
DWELL_AUT = dict(rugen=22, linum=18, diepholz=16, fochtelo=10, der=24, arjuzanx=0, gallocanta=8,
                 matsalu=18, biebrza=10, hortobagy=30, polesia=12, askania=22, anatolie=8, hula=0)
DWELL_SPR = dict(gallocanta=6, der=14, arjuzanx=0, diepholz=10, linum=8, fochtelo=6, hornborga=16,
                 hortobagy=12, biebrza=6, matsalu=8, askania=10, polesia=6, anatolie=5, hula=0)
GROUPS = [
    (0.24, 0, ["b_scan", "rugen", "diepholz", "der", "gallocanta", "extremadura"],
              ["extremadura", "gallocanta", "der", "diepholz", "hornborga", "b_scan"], 42, 205),
    (0.08, 0, ["b_scan", "rugen", "linum", "der", "arjuzanx"],
              ["arjuzanx", "der", "linum", "rugen", "hornborga", "b_scan"], 45, 200),
    (0.10, 0, ["b_ger", "linum", "der", "gallocanta", "extremadura"],
              ["extremadura", "gallocanta", "der", "linum", "b_ger"], 62, 196),
    (0.05, 0, ["b_nl", "fochtelo", "der", "arjuzanx"],
              ["arjuzanx", "der", "fochtelo", "b_nl"], 70, 205),
    (0.10, 0, ["b_balt", "rugen", "diepholz", "der", "gallocanta", "extremadura"],
              ["extremadura", "gallocanta", "der", "diepholz", "rugen", "b_balt"], 40, 200),
    (0.15, 1, ["b_fin", "matsalu", "biebrza", "hortobagy", "ichkeul"],
              ["ichkeul", "hortobagy", "biebrza", "matsalu", "b_fin"], 35, 210),
    (0.08, 1, ["b_balt", "biebrza", "hortobagy", "algerije"],
              ["algerije", "hortobagy", "biebrza", "b_balt"], 45, 205),
    (0.12, 2, ["b_rus", "polesia", "askania", "anatolie", "hula"],
              ["hula", "anatolie", "askania", "polesia", "b_rus"], 38, 215),
    (0.08, 2, ["b_rus", "polesia", "askania", "anatolie", "hula", "ethiopie"],
              ["ethiopie", "hula", "anatolie", "askania", "polesia", "b_rus"], 32, 195),
]


def haversine_km(a, b):
    lon1, lat1, lon2, lat2 = map(math.radians, (*a, *b))
    h = math.sin((lat2 - lat1) / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2
    return 6371 * 2 * math.asin(math.sqrt(h))


def chaikin(coords, n=3):
    pts = np.asarray(coords, float)
    for _ in range(n):
        q = 0.75 * pts[:-1] + 0.25 * pts[1:]
        r = 0.25 * pts[:-1] + 0.75 * pts[1:]
        pts = np.vstack([pts[:1], np.column_stack([q, r]).reshape(-1, 2), pts[-1:]])
    return pts


def site_point(key, spread=1.0):
    s = SITES[key]
    return (s["lon"] + RNG.normal(0, s["sx"] * spread), s["lat"] + RNG.normal(0, s["sy"] * spread))


def flight(a, b, t0, path, ts):
    """Vlucht a->b als gebogen boog (kwadratische bezier), 350-600 km/dag."""
    d = haversine_km(a, b)
    speed = RNG.uniform(350, 600)
    dur = max(0.4, d / speed)
    ax, ay = a; bx, by = b
    mx, my = (ax + bx) / 2, (ay + by) / 2
    nx, ny = -(by - ay), (bx - ax)
    k = RNG.normal(0, 0.08)
    cx, cy = mx + nx * k, my + ny * k
    steps = max(3, int(dur * 6))
    for i in range(1, steps + 1):
        u = i / steps
        x = (1 - u) ** 2 * ax + 2 * (1 - u) * u * cx + u ** 2 * bx
        y = (1 - u) ** 2 * ay + 2 * (1 - u) * u * cy + u ** 2 * by
        path.append((x, y)); ts.append(t0 + u * dur)
    return t0 + dur


def stay(pos, t0, t1, path, ts, jitter=0.03):
    """Verblijf: om de ~2 dagen een kleine verplaatsing (foerageren/slaapplaats)."""
    t = t0
    while t + 2 < t1:
        t += RNG.uniform(3, 5)
        path.append((pos[0] + RNG.normal(0, jitter), pos[1] + RNG.normal(0, jitter))); ts.append(t)
    path.append(pos); ts.append(t1)


def simulate_bird(group):
    _, fw, aut, spr, dep, sdep = group
    path, ts = [], []
    home = site_point(aut[0])
    t = float(np.clip(RNG.normal(dep, 9), 5, 120))
    path.append(home); ts.append(0.0)
    stay(home, 0.0, t, path, ts, 0.05)
    visits = []  # (site, t_in, t_out)
    pos = home
    for key in aut[1:]:
        nxt = site_point(key, 0.6)
        t = flight(pos, nxt, t, path, ts)
        pos = nxt
        if key == aut[-1]:
            winter_key, t_win = key, t
            break
        dwell = max(1.0, RNG.normal(DWELL_AUT.get(key, 8), DWELL_AUT.get(key, 8) * 0.35))
        visits.append((key, t, t + dwell))
        stay(pos, t, t + dwell, path, ts)
        t += dwell
    # overwintering; terugtrek
    t_spr = float(max(t + 10, RNG.normal(sdep, 8)))
    visits.append((winter_key, t_win, t_spr))
    stay(pos, t, t_spr, path, ts, 0.08)
    t = t_spr
    for key in spr[1:]:
        nxt = home if key == spr[-1] else site_point(key, 0.6)
        t = flight(pos, nxt, t, path, ts)
        pos = nxt
        if key == spr[-1]:
            break
        dwell = max(1.0, RNG.normal(DWELL_SPR.get(key, 6), DWELL_SPR.get(key, 6) * 0.35))
        visits.append((key, t, t + dwell))
        stay(pos, t, t + dwell, path, ts)
        t += dwell
    stay(pos, t, 365.0, path, ts, 0.05)
    # sorteer/clip tijden
    out_p, out_t = [], []
    for p, tt in zip(path, ts):
        if tt > 365: break
        if out_t and tt <= out_t[-1]: tt = out_t[-1] + 0.01
        out_p.append([round(p[0], 2), round(p[1], 2)]); out_t.append(round(tt, 2))
    return fw, out_p, out_t, visits


def in_flight_mask(ts, path, days):
    """Per dag: is de vogel op die dag >150 km verplaatst?"""
    arr_t = np.asarray(ts); arr_p = np.asarray(path)
    lon = np.interp(days, arr_t, arr_p[:, 0]); lat = np.interp(days, arr_t, arr_p[:, 1])
    lon2 = np.interp(days + 1, arr_t, arr_p[:, 0]); lat2 = np.interp(days + 1, arr_t, arr_p[:, 1])
    dist = np.hypot((lon2 - lon) * np.cos(np.radians(lat)), lat2 - lat) * 111
    return dist > 150


def main(src, dst):
    src, dst = Path(src), Path(dst)
    dst.mkdir(parents=True, exist_ok=True)

    # --- landen ---------------------------------------------------------
    landen = gpd.read_file(src / "ne_50m_admin_0_countries.geojson")
    landen = landen[["NAME_NL", "ISO_A2_EH", "geometry"]].rename(columns={"NAME_NL": "naam", "ISO_A2_EH": "iso"})
    landen = gpd.clip(landen, BBOX)
    landen["geometry"] = landen.geometry.simplify(0.03, preserve_topology=True)
    landen = landen[~landen.is_empty]
    landen.to_file(dst / "landen.geojson", driver="GeoJSON", COORDINATE_PRECISION=3)

    # --- water ----------------------------------------------------------
    meren = gpd.read_file(src / "ne_50m_lakes.geojson")[["name", "geometry"]]
    meren["type"] = "meer"
    rivieren = gpd.read_file(src / "ne_50m_rivers_lake_centerlines.geojson")[["name", "geometry"]]
    rivieren["type"] = "rivier"
    water = gpd.GeoDataFrame(pd_concat([meren, rivieren]), crs="EPSG:4326")
    water = gpd.clip(water, BBOX)
    water["geometry"] = water.geometry.simplify(0.02)
    water.to_file(dst / "water.geojson", driver="GeoJSON", COORDINATE_PRECISION=3)

    # --- rustplaatsen ---------------------------------------------------
    rows = []
    for k, s in SITES.items():
        rows.append(dict(id=k, naam=s["naam"], land=s["land"], rol=s["rol"], piek=s["piek"] or 0,
                         info=s["info"], geometry=Point(s["lon"], s["lat"])))
    sites = gpd.GeoDataFrame(rows, crs="EPSG:4326")
    sites.to_file(dst / "rustplaatsen.geojson", driver="GeoJSON")

    # --- flyways + corridors (afstand en buffer in EPSG:3035) -------------
    fl_rows = []
    for f in FLYWAYS:
        pts = [(SITES[k]["lon"], SITES[k]["lat"]) for k in f["hartlijn"]]
        fl_rows.append(dict(id=f["id"], naam=f["naam"], kleur=f["kleur"], geometry=LineString(chaikin(pts, 4))))
    fly = gpd.GeoDataFrame(fl_rows, crs="EPSG:4326")
    fly_laea = fly.to_crs(LAEA)
    fly["lengte_km"] = (fly_laea.length / 1000).round(0).astype(int)
    fly.to_file(dst / "flyways.geojson", driver="GeoJSON", COORDINATE_PRECISION=4)
    corr = fly_laea.copy()
    corr["geometry"] = fly_laea.buffer(150_000, quad_segs=16).simplify(5_000)
    corr = corr.to_crs("EPSG:4326")
    corr["lengte_km"] = fly["lengte_km"]
    corr.to_file(dst / "corridors.geojson", driver="GeoJSON", COORDINATE_PRECISION=3)

    # --- simulatie -------------------------------------------------------
    shares = np.array([g[0] for g in GROUPS]); shares /= shares.sum()
    gidx = RNG.choice(len(GROUPS), size=N_BIRDS, p=shares)
    birds, days = [], np.arange(366, dtype=float)
    occ = {k: np.zeros(366) for k, s in SITES.items() if s["rol"] != "broed"}
    moving = np.zeros(366)
    moving_fw = np.zeros((3, 366))
    for gi in gidx:
        fw, p, t, visits = simulate_bird(GROUPS[gi])
        birds.append(dict(f=int(fw), p=p, t=t))
        for key, tin, tout in visits:
            if key in occ:
                occ[key][(days >= tin) & (days < tout)] += 1
        m = in_flight_mask(t, p, days)
        moving += m; moving_fw[fw] += m
    # schaal bezetting per site naar piektelling
    occ_out = {}
    for k, arr in occ.items():
        mx = arr.max()
        occ_out[k] = [int(round(v / mx * SITES[k]["piek"], -2)) if mx else 0 for v in arr]
    json.dump(dict(start="08-01", eenheid="dag", n=N_BIRDS,
                   flyways=[dict(id=f["id"], naam=f["naam"], kleur=f["kleur"],
                                 lengte_km=int(fly.loc[f["id"], "lengte_km"])) for f in FLYWAYS],
                   birds=birds),
              open(dst / "tracks.json", "w"), separators=(",", ":"))
    json.dump(dict(bezetting=occ_out,
                   in_trek=[round(v / N_BIRDS, 3) for v in moving],
                   in_trek_flyway=[[round(v / N_BIRDS, 3) for v in r] for r in moving_fw]),
              open(dst / "bezetting.json", "w"), separators=(",", ":"))

    print(fly[["naam", "lengte_km"]].to_string(index=False))
    print("vogels:", N_BIRDS, "punten:", sum(len(b["p"]) for b in birds))


def pd_concat(frames):
    import pandas as pd
    return pd.concat(frames, ignore_index=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
