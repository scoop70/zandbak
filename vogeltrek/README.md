# Vogeltrek – kraanvogel (Grus grus)

`../vogeltrek.html` toont een geanimeerde kaart van de kraanvogeltrek door Europa (MapLibre GL + deck.gl).

| Bestand | Inhoud |
|---|---|
| `build_vogeltrek.py` | GeoPandas-pipeline (Natural Earth → landen/water; flyways, corridors in EPSG:3035; simulatie) |
| `landen.geojson`, `water.geojson` | Natural Earth 1:50m, geclipt en vereenvoudigd |
| `rustplaatsen.geojson` | broed-, rust- en overwinteringsgebieden met orde-van-grootte piektellingen |
| `flyways.geojson`, `corridors.geojson` | hartlijnen (Chaikin) + 150 km buffer |
| `tracks.json` | 640 **gemodelleerde** individuen (geen GPS-data), tijd in dagen sinds 1 aug |
| `bezetting.json` | gemodelleerde bezetting per rustplaats per dag, aandeel in trek per flyway |

Live in de browser: GBIF Occurrence API (alle bronnen of alleen eBird Observation Dataset),
GBIF hex-dichtheidstegels, NOAA NCEI DEM Global Mosaic hillshade (ETOPO 2022), terugval Esri World Hillshade Dark.

Opnieuw bouwen:

```bash
python build_vogeltrek.py <map-met-natural-earth-geojson> .
```
